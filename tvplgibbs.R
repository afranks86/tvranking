############################################################################
## [lambda_st, a_st, rho_st, stats] = tvplack(X, K, a, rho, rho_param,
## N_Gibbs, N_burn)
##
## Gibbs sampler for dynamic Plackett-Luce model
##
##########
## INPUTS
##
## X is a cell of size N where each entry is a Nx3 matrix, where each row contains
##   Column 1:  individual ID (1 through K)
##   Column 2:  contest ID (1 through n)
##   Column 3:  rank 
## K is the total number of individuals to be ranked (may be higher than the
## number of elements in X)
## a is the shape parameter for the gamma prior (default a = 1)
##  If a<0, then it is estimated with a vague prior
## rho is the parameter tuning correlations between skill ratings at
## different times
## rho_param is a matrix with parameters of the gamma prior for different
## values of rho
## N_Gibbs: Number of Gibbs iterations
## N_burn: Number of burn-in iterations
##########
## OUTPUTS
## lambda_st gives the values of the skills at each iteration
## a_st gives the values of the shape parameter at each iteration
## stats is a structure with some summary statistics on the parameters
###########################################################################
## Reference:
## L. Bornn and F. Caron. Dynamic generalized Bradley Terry models, 2011.
##
## August 2011
## Author: F. Caron (INRIA Bordeaux Sud-Ouest)
## Francois.Caron@inria.fr
## http://www.math.u-bordeaux1.fr/~fcaron/
###########################################################################

tvplgibbs <- function(ranks, bmask.list=NULL,
                          kappa=1, nu=1,omega=Inf,
                          rho, rho_param, Ngibbs, Nburn=1,
                          init.list=NULL,sample.skip=10,drawRho=FALSE){


    ## Number of time steps, T, and number of individuals K
    T <- ncol(ranks)
    K <- nrow(ranks)

    ## init rho_mask for rho estimation
    ## rho is a K x T matrix (how many different values of rho are there?)
    rho_val <- unique(as.vector(rho))
    estimate_rho <- rep(1,length(rho_val))
    rho_mask <- array(0,dim=c(K,T-1,length(rho_val)))
    for(i in 1:length(rho_val)){
        rho_mask[,,i] <- rho==rho_val[i]
    }

    if(a<0){
        estimate_a <- TRUE
        a <- -a
    }
    else{
        estimate_a <- FALSE
    }

    if(is.null(bmask.list)){
        bmask.list[["ALL"]] <- matrix(TRUE,nrow=K,ncol=T,dimnames=list(rownames(wins[[1]]),1:T))
    }
    
        
    ## Initialization
    if(is.null(init.list)){
        init.list <- init_bthometies(wins,losses,ties,a,bmask.list)
    }
    lambda <- init.list$lambda
    indices1 <- init.list$indices1
    indices2 <- init.list$indices2
    theta <- init.list$theta
    alpha <- init.list$alpha
    at <- init.list$at
    bi <- init.list$bi
    bmat <- init.list$bmat
    b.groups <- matrix(0,nrow=length(bmask.list),ncol=T,dimnames=list(names(bmask.list),1:T))


    ## Store every skip-th value
    Nsamples <- (Ngibbs-Nburn)/sample.skip
    a.samps <- rep(1,Nsamples)
    b.samps <- array(0,dim=c(length(bmask.list),T,Nsamples),dimnames=list(names(bmask.list),1:T,1:Nsamples))
    rho.samps <- matrix(0,nrow=length(rho_val),ncol=Nsamples)
    lambda.samps <- array(0,dim=c(K,T,Nsamples))
    

    lambda.samps[,,1] <- matrix(0,nrow=K,ncol=T)
    for(t in 1:T){
        lambda.samps[,t,1] <- lambda[,t]/sum(lambda[,t])*K 
    }

    theta.samps[1] <- init.list$theta
    alpha.samps[1] <- init.list$alpha
    a.samps[1] <- a
    b.samps[,,1] <- bmat

    ## Init latent variables, V

    V <- matrix(rpois(K*(T-1),as.vector(bi[,1:(T-1)]*rho*lambda[,1:(T-1)])),nrow=K,ncol=(T-1))

    ## Init latent variables W, to mean
    if(omega!=Inf) {
        W <- matrix(round(as.vector(nu*omega*bmat[,1:(T-1)])),nrow=length(bmask.list),ncol=(T-1),dimnames=list(names(bmask.list),1:(T-1)))
    } else {
        W <- NULL
    }
    
    ## Run sampler for NGibbs iterations
    for(i in 1:Ngibbs){

        if(i%%1==0)
            print(sprintf("i = %i",i))

        ## sample total bi
        B <- rep(1,T)
        B[1] <- sum(rgamma(nrow(bmat),kappa,nu))
        bmat[,1] <- B[1]*bmat[,1]/sum(bmat[,1])
        if(omega==Inf) {
            bmat[,2:T] <- bmat[,1]
        } else {
            Wtmp <- rpois(nrow(bmat),nu*omega*bmat[,1])
            for(t in 2:T) {
                B[t] <- sum(rgamma(nrow(bmat),kappa+Wtmp,nu*(1+omega)))
                bmat[,t] <- B[t]*bmat[,t]/sum(bmat[,t])
                if(t<T)
                    Wtmp <- rpois(nrow(bmat),nu*omega*bmat[,t])
            }
        }
        for(grp in names(bmask.list)) {
            bmask <- bmask.list[[grp]]
            bi[bmask] <- t(t(bmask)*bmat[grp,])[bmask]
        }


        
        ## Sample total weights at each time from prior
        S <- rep(1,T)

        lambda.prior <- rgamma(K,a,bi[,1])
        S[1] <- sum(lambda.prior)
        lambda[,1] <- S[1]*lambda[,1]/sum(lambda[,1])
        Vtmp <- rpois(K,bi[,2]*rho[,1]*lambda.prior)
        for(t in 2:T){
            lambda.prior <- rgamma(K,a+Vtmp,bi[,t]*(1+rho[,t-1]))
            S[t]  <- sum(lambda.prior)
            lambda[,t] <- lambda[,t]*S[t]/sum(lambda[,t])
            if(t<T)
                Vtmp <- rpois(K,bi[,t]*rho[,(t-1)]*lambda.prior)
        }

        ## Sample the latent V
        V <- sampleV(V,lambda,rho,a,bi)

        ## Update bthometies
        ## update <- update_bthometies(wins,losses,ties,indices1,indices2,T,K,V,lambda,alpha,theta,at,bi)
        ## lambda <- update$lambda
        ## alpha <- update$alpha
        ## theta <- update$theta

        ## Sample rho given lambda, a
        for(k in 1:length(rho_val)){
            if(estimate_rho[k]&drawRho){
                rho.draw <- sampleRho(lambda,rho,rho_val[k],rho_mask[,,k],a,bi,rho_param[k,1],rho_param[k,2])
                rho <- rho.draw$rho
                rho_val[k] <- rho.draw$rho_val
            }
        }

        ## sample bi

        bout <- samplebi(bi,bmask.list,bmat,lambda=lambda,rho=rho,V=V,W=W,T=T,kappa=kappa,nu=nu,omega=omega)
        bi <- bout$bi
        bmat <- bout$bmat
        W <- bout$W
        
        #print(rho[1,])
        ## sample a given rho, lambda
        if(estimate_a)
            a <- samplea(rho,lambda,a,b)

        ## Store current sample if not burn-in
        if(i>Nburn){
            ind <- ceiling((i-Nburn)/sample.skip)
            lambda.samps[,,ind] <- lambda
            theta.samps[ind] <- theta
            alpha.samps[ind] <- alpha
            a.samps[ind] <- a
            rho.samps[,ind] <- rho_val
            b.samps[,,ind] <- bmat
        }
        
    }

    list(lambda.samps=lambda.samps,theta.samps=theta.samps,alpha.samps=alpha.samps,
         a.samps=a.samps,rho.samps=rho.samps,b.samps=b.samps)
    
}

sampleV <- function(V,lambda,rho,a,bi,N.MH=5){

    K <- nrow(lambda)
    T <- ncol(lambda)

    for(nn in 1:N.MH){  ##number of metropolis proposals

        Vnew <- rpois(K*(T-1),lambda=as.vector(bi[,2:T]*rho*lambda[, 1:(T-1)]))

        ## Need this line?
        Vnew <- matrix(Vnew,nrow=K,ncol=(T-1))

        ## Double check this!!!
        u <- runif(K*(T-1))
        logaccept <- (-V+Vnew)*log(bi[,2:T]*(1+rho)*lambda[, 2:T]) + lgamma(a+V) - lgamma(a+Vnew)
        
        accept <- as.vector(exp(logaccept))
                                   
        V[u<accept] <- Vnew[u<accept]
    }

    ## Check this too!
    matrix(V,nrow=K,ncol=(T-1))
}

samplea <- function(rho,lambda,a,b,T,N.MH=2){

    for(iter in 1:N.MH){
        print("TODO")
        return(NULL)
    }
    
}

## For now assume constant bi
## b_1,i ~ gamma(kappa,nu)
## b_t,i|W_t-1 ~ Gam(kappa+W_t-1,nu*(1+w_t-1))
samplebi <- function(bi,bmask.list,bmat,lambda,rho,V,W,T,kappa,nu,omega=Inf,N.MH=5){

    ## For non-time-varying bi
    if(omega==Inf){
        for(grp in names(bmask.list)){

            bmask.full <- bmask.list[[grp]]
            bmask.first <- bmask.full[,1]
            bmask <- bmask.full[,2:T]
            
            bcur.rho <- rho[bmask]
            cur.lambda1 <- (lambda[,1:(T-1)])[bmask]
            cur.lambda2 <- (lambda[,2:T])[bmask]
            
            b.alpha <- sum(2*V[bmask]+a)+sum(bmask.first)*a+kappa
            b.beta <- sum(bcur.rho*cur.lambda1+(1+bcur.rho)*cur.lambda2)+sum(lambda[bmask.first,1])+nu

            b.new <- rgamma(1,b.alpha,b.beta)
            bmat[grp,] <- b.new
            bi[bmask.full] <- b.new
        }
    } else {

        for(grp in names(bmask.list)){

            bmask <- bmask.list[[grp]]

            ## Sample W
            Wvec.cur <- W[grp,]
            for(nn in 1:N.MH){

                ## independence proposal
                Wvec.new <- rpois((T-1),lambda=as.vector(nu*omega*bmat[grp,1:(T-1)]))

                ## Double check this!!!
                u <- runif(T-1)
                logaccept <- (Wvec.new-Wvec.cur)*
                    (log(bmat[grp,1:(T-1)]*nu*omega) +
                     log(nu*(1+omega)) +
                     log(bmat[grp,2:T])) +
                    lgamma(kappa+Wvec.cur) - lgamma(kappa+Wvec.new) + lgamma(Wvec.cur+1) - lgamma(Wvec.new+1)
        
                accept <- as.vector(exp(logaccept))
                                   
                Wvec.cur[u<accept] <- Wvec.new[u<accept]
            }
            
            W[grp,] <- Wvec.cur

            b.new <- numeric(T)
            ## T=1
            b.new[1] <- rgamma(1,a*sum(bmask[,1])+Wvec.cur[1]+kappa,
                               sum(lambda[bmask[,1],1])+nu*(1+omega))

            ## T = T
            b.new[T] <- rgamma(1,sum(2*V[bmask[,T],T-1]+a)+kappa+Wvec.cur[T-1],
                               sum(lambda[bmask[,T],T]*(1+rho[bmask[,T],T-1])+
                                   lambda[bmask[,T],T-1]*rho[bmask[,T],T-1])+nu*(1+omega))


            ## 2:(T-1)
            b.new[2:(T-1)] <- rgamma(T-2, colSums((2*V[,1:(T-2)]+a)*bmask[,2:(T-1)])+kappa+Wvec.cur[1:(T-2)]+Wvec.cur[2:(T-1)], colSums((rho[,1:(T-2)]*lambda[,1:(T-2)]+lambda[,2:(T-1)]*(1+rho[,1:(T-2)]))*bmask[,2:(T-1)])+nu*(1+2*omega))

            bmat[grp,] <- b.new
            bi[bmask] <- t(t(bmask)*b.new)[bmask]

        }
    }

    list(bi=bi,bmat=bmat,W=W)
}

sampleRho <- function(lambda,rho,rho_val,rho_mask,a,bi,rho_a, rho_b,N.MH=2){

    T <- ncol(lambda)
    for(iter in 1:N.MH){ ## Number of MH iterations

        ## propose
        rho_val_new <- rho_val*exp(.1*rnorm(1))
        rho_new <- rho
        rho_new[as.logical(rho_mask[])] <- rho_val_new

        temp <- 2*sqrt(bi[,2:T]*(1+rho)*bi[,2:T]*rho*lambda[, 1:(T-1)]*lambda[, 2:T])
        temp_new <- 2*sqrt(bi[,2:T]*(1+rho_new)*bi[,2:T]*rho_new*lambda[, 1:(T-1)]*lambda[, 2:T])
        
        lograte <- (rho_new - rho)*(lambda[, 1:(T-1)] + lambda[, 2:T])+(a+1)/2*log(bi[,2:T]*rho)-(a-1)/2*log(rho)-(a+1)/2*log(bi[,2:T]*rho_new)+(a-1)/2*log(rho_new)+log(besselI(as.numeric(temp),a-1,1))-log(besselI(as.numeric(temp_new),a-1,1)) + temp - temp_new  
        ## correction for using besseli(:,:,1)
        
        lograte <- -sum(sum(lograte))+log(rho_val_new) - log(rho_val)+dgamma(rho_val_new, rho_a, rho_b,log=TRUE)-dgamma(rho_val, rho_a, rho_b,log=TRUE)
        
        if(runif(1)<exp(lograte)){ ## If accept
            rho <- rho_new
            rho_val <- rho_val_new
        }
    }

    list(rho=rho,rho_val=rho_val)
}

update_pl <- function(ranks,T,K,V,lambda,at,bi) {

    ## Sample Z given lambda
    Z <- vector("list",T)
    bdat <- matrix(0,nrow=K,ncol=T)

    for(t in 1:T){
        indices <- order(which(!is.na(lambda[,t]))) ## preprocess this?
        Z[[t]] <- rexp(cumsum(lambda[indices,1])) ## todo order by ranks
    }

    tmp <- sapply(Z,cumsum)

    ## Sample lambda given Z,U,theta,alpha
    lambda[,1] <- rgamma(K,at[,1]+V[,1],bi[,1](1+rho[,1]))

    lambda[,2:(T-1)] <- rgamma(K*(T-2),at[,2:(T-1)]+V[,2:(T-1)]+V[,1:(T-2)],
                               bdat[,2:(T-1)]+bi[,2:(T-1)]*rho[,2:(T-1)]+bi[,1:(T-2)]*(1+rho[,1:(T-2)]))

    lambda[,T] <- rgamma(K,at[,T]+V[,T-1],bdat[,T]+bi[,T-1]*(1+rho[,T-1]))
        
    
    lambda
}

