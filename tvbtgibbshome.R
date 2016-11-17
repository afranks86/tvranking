# NEW MODEL WITH POISSON LATENT VARIABLES


###########################################################################
# function [pi_st, theta_st, alpha_st, a_st, stats] = tvbtgibbshometies(w, l,
# t, a, N_Gibbs, N_burn)
#
# Gibbs sampler for Bradley-Terry models with home advantage and ties
#
# Requires the statistics toolbox
##########
# INPUTS
#
# wins[[k]][i, j] is the number of times i beats j when i is at home at time k
# losses[[k]][i, j] is the number of times i looses to j when i is at home at time k
# ties[[k]][i, j] is the number of times i ties j when i is at home at time k
# a is the shape parameter for the gamma prior (default a = 1)
#  If a<0, then it is estimated with a vague prior
# rho is the parameter tuning correlations between skill ratings at
# different times
# rho_param is a length(rho_val)x2 matrix with parameters of the gamma prior for different
# values of rho
# values of rho
# N_Gibbs: Number of Gibbs iterations
# N_burn: Number of burn-in iterations
##########
# OUTPUTS
# lambda_st gives the values of the skills at each iteration
# theta_st gives the values of the ties parameter at each iteration
# alpha_st gives the values of the home advantage parameter at each iteration
# a_st gives the values of the shape parameter at each iteration
# stats is a structure with some summary statistics on the parameters
###########################################################################
# Reference:
# L. Bornn and F. Caron. Dynamic generalized Bradley Terry models, 2011.
#
# April 2011
# Author: Franks, Alexander and Caron (INRIA Bordeaux Sud-Ouest), Alex Franks
# Francois.Caron@inria.fr
# http://www.math.u-bordeaux1.fr/~fcaron/
###########################################################################

tvbtgibbshome <- function(wins, losses, ties, a, bmask.list=NULL,
                          kappa=1, nu=1,omega=Inf,
                          rho, rho_param, Ngibbs, Nburn=1,
                          init.list=NULL,sample.skip=10, draw=c()){



    ## whcich parameters to sample (useful for debugging)
    draw <- c(draw, V=TRUE, b=TRUE, rho=FALSE, a=FALSE)
    draw <- draw[unique(names(draw))]
        
    ## Number of time steps, T, and number of individuals K
    T <- length(wins)
    K <- nrow(wins[[1]])

    ## init rho_mask for rho estimation
    ## rho is a K x T matrix (how many different values of rho are there?)
    rho_val <- unique(as.vector(rho))
    estimate_rho <- rep(1,length(rho_val))
    rho_mask <- array(0,dim=c(K,T-1,length(rho_val)))
    for(i in 1:length(rho_val)){
        rho_mask[, , i] <- rho==rho_val[i]
    }

    if(is.null(bmask.list)){
        bmask.list[["ALL"]] <- matrix(TRUE,nrow=K,ncol=T,dimnames=list(rownames(wins[[1]]),1:T))
    }
    
        
    ## Initialization
    if(is.null(init.list)){
        init.list <- init_bthometies(wins, losses, ties, a, bmask.list)
    }
    
    lambda <- init.list$lambda
    theta <- init.list$theta
    alpha <- init.list$alpha
    wltStat <- init.list$wltStat

    indices1 <- lapply(1:T, function(t) summary(wins[[t]] + ties[[t]]))
    indices2 <- lapply(1:T, function(t) summary(losses[[t]]+ties[[t]]))

    bi <- init.list$bi
    bmat <- init.list$bmat
    
    ## Store every skip-th value
    Nsamples <- (Ngibbs-Nburn)/sample.skip
    theta.samps <- rep(0,Nsamples)
    alpha.samps <- rep(0,Nsamples)
    a.samps <- rep(1,Nsamples)
    b.samps <- array(0, dim=c(length(bmask.list), T, Nsamples),
                     dimnames=list(names(bmask.list), 1:T, 1:Nsamples))
    rho.samps <- matrix(0,nrow=length(rho_val),ncol=Nsamples)
    lambda.samps <- array(0,dim=c(K,T,Nsamples))
    rownames(lambda.samps) <- rownames(wins[[1]])
    

    lambda.samps[ , ,1] <- matrix(0, nrow=K, ncol=T)
    for(t in 1:T){
        lambda.samps[, t, 1] <- lambda[, t]
    }

    theta.samps[1] <- init.list$theta
    alpha.samps[1] <- init.list$alpha
    a.samps[1] <- a
    b.samps[, , 1] <- bmat

    ## Init latent variables, V
    if(is.null(init.list$V))
        V <- matrix(rpois(K*(T-1), as.vector(bi[, 1:(T-1)]*rho*lambda[, 1:(T-1)])),
                    nrow=K, ncol=(T-1))
    else
        V <- init.list$V

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

        for(grp in names(bmask.list)) {
            bmask <- bmask.list[[grp]]
            bi[bmask] <- t(t(bmask) * bmat[grp, ])[bmask]
        }

        S <- rgamma(1, K*a, bi[1])
        lambda <- lambda/mean(colSums(lambda))*S 

        ## Sample the latent V
        if(draw["V"])
            V <- sampleV(V, lambda, rho, a, bi)

        ## Update bthometies`
        if(draw["lambda"]) {
            update <- update_bthometies(wins, losses, ties,
                                        indices1, indices2,
                                        T, K, V, lambda, alpha, theta, wltStat, bi)
            lambda <- update$lambda
            alpha <- update$alpha
            theta <- update$theta
        }
        
        ## Sample rho given lambda, a
        if(draw["rho"]) {
            for(k in 1:length(rho_val)) {
                rho.draw <- sampleRho(lambda, rho, rho_val[k],
                                      rho_mask[, , k], a, bi,
                                      rho_param[k, 1], rho_param[k, 2])
                    rho <- rho.draw$rho
                    rho_val[k] <- rho.draw$rho_val
            }
        }

        ## sample bi
        if(draw["b"]) {
            
            bout <- samplebi(bi, bmask.list, bmat, lambda=lambda, rho=rho,
                             V=V, W=W, T=T, kappa=kappa, nu=nu, omega=omega)
            bi <- bout$bi
            bmat <- bout$bmat

            bi <- t(t(bi)/colMeans(bi)) * kappa/nu
            bmat <- t(t(bmat)/colMeans(bi)) * kappa/nu
            W <- bout$W
        }
        
        if(draw["a"])
            a <- samplea(rho,lambda,a,b)

        ## temporary
        print(mean(lambda))
        print(median(as.numeric(V)))
        
        ## Store current sample if not burn-in
        if(i > Nburn){

            ind <- ceiling((i-Nburn)/sample.skip)
            lambda.samps[,,ind] <- lambda
            theta.samps[ind] <- theta
            alpha.samps[ind] <- alpha
            a.samps[ind] <- a
            rho.samps[,ind] <- rho_val
            b.samps[,,ind] <- bmat
        }
        ## if(i %%  100 == 0)
            ## browser()
        
    }

    list(lambda.samps=lambda.samps,
         theta.samps=theta.samps,alpha.samps=alpha.samps,
         a.samps=a.samps, rho.samps=rho.samps, b.samps=b.samps)
    
}

sampleV <- function(V, lambda, rho, a, bi, N.MH=5){

    K <- nrow(lambda)
    T <- ncol(lambda)

    for(nn in 1:N.MH){  ##number of metropolis proposals

        Vnew <- rpois(K*(T-1), lambda=as.vector(bi[, 1:(T-1)]*rho*lambda[, 1:(T-1)]))
        u <- runif(K*(T-1))

        logaccept <- (Vnew - V)*log(bi[, 2:T]*(1+rho)*lambda[, 2:T]) + lgamma(a+V) - lgamma(a+Vnew)

        accept <- as.vector(exp(logaccept))

        V[u < accept] <- Vnew[u < accept]
    }

    ## Check this too!
    matrix(V, nrow=K, ncol=(T-1))
}

samplea <- function(rho, lambda, a, b, T, N.MH=2){

    for(iter in 1:N.MH){
        print("TODO")
        return(NULL)
    }
    
}

## For now assume constant bi
## b_1,i ~ gamma(kappa,nu)
## b_t,i|W_t-1 ~ Gam(kappa+W_t-1,nu*(1+w_t-1))
samplebi <- function(bi, bmask.list, bmat, lambda, rho, V, W, T,
                     kappa, nu, omega=Inf, N.MH=5){
    ## For non-time-varying bi
    if(omega==Inf){
        for(grp in names(bmask.list)){
            bmask.full <- bmask.list[[grp]]
            bmask.first <- bmask.full[, 1]
            bmask.last <- bmask.full[, T]
            bmaskMid <- bmask.full[, 2:(T-1)]
            
            b.alpha <-
                sum(V[, 1:(T-2)][bmaskMid] + V[, 2:(T-1)][bmaskMid]) +
                sum(V[bmask.first, 1] ) +
                sum(V[bmask.last, T-1]) +
                sum(bmask.full)*a
                kappa
            b.beta <- sum((1 + rho[, 1:(T-2)][bmaskMid] +
                           rho[, 2:(T-1)][bmaskMid])*lambda[, 2:(T-1)][bmaskMid]) +
                sum((1+rho[bmask.first, 1])*lambda[bmask.first, 1]) +
                sum((1+rho[bmask.last, T-1])*lambda[bmask.last, T]) + nu
                
            b.new <- rgamma(1, b.alpha, b.beta)
            bmat[grp, ] <- b.new
            bi[bmask.full] <- b.new

        }

    } else {

        for(grp in names(bmask.list)){

            bmask <- bmask.list[[grp]]

            ## Sample W
            Wvec.cur <- W[grp,]
            for(nn in 1:N.MH){

                ## independence proposal
                Wvec.new <- rpois((T-1), lambda=as.vector(nu*omega*bmat[grp,1:(T-1)]))

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

update_bthometies <- function(wins, losses, ties, indices1, indices2,
                              T, K, V, lambda, alpha, theta, wltStat, bi){

    ## Sample Z given lambda
    Z1 <- Z2 <- vector("list",T)
    bdat <- matrix(0, nrow=K, ncol=T)

    total_ties <- sum(sapply(ties,sum))
    total_wins_and_ties <- total_ties + sum(sapply(wins,sum))
    
    for(t in 1:T){
    
        N1 <- nrow(indices1[[t]])
        i1 <- indices1[[t]]$i
        j1 <- indices1[[t]]$j
        val1 <- indices1[[t]]$x

        N2 <- nrow(indices2[[t]])
        i2 <- indices2[[t]]$i
        j2 <- indices2[[t]]$j
        val2 <- indices2[[t]]$x
        
        Z1[[t]] <- sparseMatrix(i=i1, j=j1, x=rgamma(N1, shape=val1, rate=(alpha*lambda[i1, t]+theta*lambda[j1,t])), dims=c(K, K),check=FALSE)
        Z2[[t]] <- sparseMatrix(i=i2, j=j2, x=rgamma(N2, shape=val2, rate=(alpha*lambda[i2, t]+theta*lambda[j2,t])), dims=c(K, K),check=FALSE)
         
        bdat[, t] <- rowSums(Z1[[t]]) + theta*colSums(Z1[[t]]) + alpha*theta*rowSums(Z2[[t]]) + colSums(Z2[[t]])

    }

    ## Sample lambda given Z,U,theta,alpha
    if(T==1){
        lambda[, 1] <- rgamma(K, a + wltStat[, 1], bi[, 1] + bdat[, 1])
    }
    else{
        lambda[, 1] <- rgamma(K, a + wltStat[, 1] + V[, 1],
                              bdat[, 1] + bi[, 1]*(1 +rho[, 1]))
        lambda[, T] <- rgamma(K, a + wltStat[, T] + V[, T-1], bdat[, T] + bi[, T]*(1+rho[,T-1]))
        lambda[, 2:(T-1)] <-
            rgamma(K*(T-2), a + wltStat[, 2:(T-1)] + V[, 1:(T-2)] + V[, 2:(T-1)],
                   bdat[, 2:(T-1)] + bi[, 2:(T-1)]*rho[, 2:(T-1)] + bi[, 1:(T-2)]*(1+rho[, 1:(T-2)]))

    }

    ## ## sample a given blah
    ## adens <- function(aval) { 
    ##     part1 <- sum(log(aval + V[, 1])*(bdat[, 1] + bi[, 1])
    ##                  - lgamma(aval + V[, 1]) + (aval+V[, 1] - 1))
    ##     part2 <- sum(log(aval + V[, T-1])*(bdat[, T] + bi[, T]*(1+rho[, T-1]))
    ##                  - lgamma(aval + V[, T-1]) + (aval+V[, T-1] - 1))


    ##     astar <- aval + V[, 2:(T-1)] + V[, 1:(T-2)]
    ##     bstar <- bdat[,2:(T-1)] + bi[,3:T]*rho[,2:(T-1)] +
    ##         bi[,2:(T-1)]*(1+rho[,1:(T-2)])
        
    ##     part3 <- sum(log(astar)*bstar - lgamma(astar) + (astar-1))

    ##     part1 + part2 + part3

    ## }
    ## a <- at[1, 1]
    ## a_new <- rlnorm(1, log(a), .1)
    ## lograte <-  adens(a_new) - adens(a) -
    ##     dlnorm(a_new, log(a), 0.1) +
    ##     dlnorm(a, log(a_new), 0.1)
    ## if(runif(1) < exp(lograte)) # If accept
    ##     a <- a_new
    ## at[] <- a
    
    ## sample alpha given Z,U,lambda
    b_alpha <- 0
    for(t in 1:T){
        b_alpha <- b_alpha + lambda[, t]%*%(rowSums(Z1[[t]]) + theta*rowSums(Z2[[t]]))
    }

    alpha <- rgamma(1, total_wins_and_ties + 20, b_alpha + 10)
    
    
    ## sample theta given Z, U,lambda with MH step
    theta_new <- rnorm(1, theta, .1) # proposal ## use truncated normal to improve eff
    if( theta_new > 1){
        temp <- 0
        for( t in 1:T)
            temp <- temp + lambda[, t] %*% (colSums(Z1[[t]]) + alpha*colSums(Z2[[t]]))

        lograte <-  total_ties*log(theta_new^2-1) - theta_new*temp-total_ties*log(theta^2 - 1+.Machine$double.eps) + theta*temp
        if(runif(1) < exp(lograte)) # If accept
            theta <- theta_new
    }

    list(lambda=lambda, alpha=alpha, theta=theta)
}

## Visualize lambdas over time
plotLambdas <- function(lambda.mean,teams){

    matplot(t(lambda.mean),type="l",xlim=c(0,ncol(lambda.mean)+20),lwd=2)
    legend("topright",teams,lty=rep(1:5,length.out=length(teams)),col=rep(1:6,length.out=length(teams)),cex=.7,lwd=2)
    grid(nx=20)
}
