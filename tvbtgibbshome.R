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
# Author: F. Caron (INRIA Bordeaux Sud-Ouest)
# Francois.Caron@inria.fr
# http://www.math.u-bordeaux1.fr/~fcaron/
###########################################################################

tvbtgibbshome <- function(wins, losses, ties, a, rho, rho_param, Ngibbs, Nburn=1,
                          init.list=NULL,sample.skip=10,drawRho=FALSE){

    
    ## Number of time steps, T, and number of individuals K
    T <- length(wins)
    K <- nrow(wins[[1]])

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

    ## How to set b?
    b <- 2
    
    ## Initialization
    if(is.null(init.list)){
        init.list <- init_bthometies(wins,losses,ties,a)
    }
    lambda <- init.list$lambda
    indices1 <- init.list$indices1
    indices2 <- init.list$indices2
    theta <- init.list$theta
    alpha <- init.list$alpha
    at <- init.list$at
    
    ## Store every skip-th value
    Nsamples <- (Ngibbs-Nburn)/sample.skip
    theta.samps <- rep(0,Nsamples)
    alpha.samps <- rep(0,Nsamples)
    a.samps <- rep(1,Nsamples)
    rho.samps <- matrix(0,nrow=length(rho_val),ncol=Nsamples)
    lambda.samps <- array(0,dim=c(K,T,Nsamples))

    lambda.samps[,,1] <- matrix(0,nrow=K,ncol=T)
    for(t in 1:T){
        lambda.samps[,t,1] <- lambda[,t]/sum(lambda[,t]) 
    }

    theta.samps[1] <- init.list$theta
    alpha.samps[1] <- init.list$alpha
    a.samps[1] <- a

    ## Init latent variables, V

    V <- matrix(rpois(K*(T-1),as.vector(rho*lambda[,1:(T-1)])),nrow=K,ncol=(T-1))

    ## Run sampler for NGibbs iterations
    for(i in 1:Ngibbs){

        if(i%%1==0)
            print(sprintf("i = %i",i))

        ## Sample total weights at each time
        S <- rep(1,T)
        S[1] <- rgamma(1,K*a,b)
        lambda[,1] <- S[1]*lambda[,1]/sum(lambda[,1])
        for(t in 2:T){
            ## Use a normal approximation to the poisson distribution to speed up code
            temp <- rnorm(1,rho[1, t-1]*S[t-1], sqrt(rho[1, t-1]*S[t-1]))
            temp[temp<0] <- 0
            S[t]  <- rgamma(1,K*a + temp, (b + rho[1, t-1]))
            lambda[,t] <- lambda[,t]*S[t]/sum(lambda[,t])
        }


        ## Sample the latent V
        V <- sampleV(V,lambda,rho,a,b)
        
        ## TODO: Update bthometies
        update <- update_bthometies(wins,losses,ties,indices1,indices2,T,K,V,lambda,alpha,theta,at,b)
        lambda <- update$lambda
        alpha <- update$alpha
        theta <- update$theta

        ## Sample rho given lambda, a
        for(k in 1:length(rho_val)){
            if(estimate_rho[k]&drawRho){
                rho.draw <- sampleRho(lambda,rho,rho_val[k],rho_mask[,,k],a,b,rho_param[k,1],rho_param[k,2])
                rho <- rho.draw$rho
                rho_val[k] <- rho.draw$rho_val
            }
        }
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
        }
        
    }

    list(lambda.samps=lambda.samps,theta.samps=theta.samps,alpha.samps=alpha.samps,
         a.samps=a.samps,rho.samps=rho.samps)
    
}

sampleV <- function(V,lambda,rho,a,b,N.MH=5){

    K <- nrow(lambda)
    T <- ncol(lambda)

    for(nn in 1:N.MH){  ##number of metropolis proposals
        
        ## Use normal approximation to the poisson
        ## Vnew <- poissrnd(rho.*(lambda(:, 1:N-1)))

        Vnew <- rnorm(K*(T-1),as.vector(rho*lambda[, 1:(T-1)]), as.vector(sqrt(rho*lambda[, 1:(T-1)])))
        Vnew[Vnew<0] <- 0

        ## Need this line?
        Vnew <- matrix(Vnew,nrow=K,ncol=(T-1))

        ## Double check this!!!
        u <- runif(K*(T-1))
        logaccept <- (-V+Vnew)*log((b+rho)*lambda[, 2:T]) + lgamma(a+V) - lgamma(a+Vnew)
        
        accept <- as.vector(exp(logaccept))

        V[u<accept] <- Vnew[u<accept]
    }

    ## Check this too!
    matrix(V,nrow=K,ncol=(T-1))
}

samplea <- function(rho,lambda,a,b,T,N.MH=2){

    for(iter in 1:N.MH){

        m_a <- 100
        var_a <- 100^2

        ## propose new from lognormal
        a_new <- a*exp(.1*rnorm())

        temp <- 2*sqrt((b+rho)*rho*lambda[,1:(T-1)]*lambda[,2:T])

        if(runif(1)<exp(lograte)){
            a <- a_new
        }
    }
    
}

sampleRho <- function(lambda,rho,rho_val,rho_mask,a,b,rho_a, rho_b,N.MH=2){

    T <- ncol(lambda)
    for(iter in 1:N.MH){ ## Number of MH iterations

        ## propose
        rho_val_new <- rho_val*exp(.1*rnorm(1))
        rho_new <- rho
        rho_new[as.logical(rho_mask[])] <- rho_val_new

        temp <- 2*sqrt((b+rho)*rho*lambda[, 1:(T-1)]*lambda[, 2:T])
        temp_new <- 2*sqrt((b+rho_new)*rho_new*lambda[, 1:(T-1)]*lambda[, 2:T])
        
        lograte <- (rho_new - rho)*(lambda[, 1:(T-1)] + lambda[, 2:T])+(a+1)/2*log(b+rho)-(a-1)/2*log(rho)-(a+1)/2*log(b+rho_new)+(a-1)/2*log(rho_new)+log(besselI(as.numeric(temp),a-1,1))-log(besselI(as.numeric(temp_new),a-1,1)) + temp - temp_new  
        ## correction for using besseli(:,:,1)
        
        lograte <- -sum(sum(lograte))+log(rho_val_new) - log(rho_val)+dgamma(rho_val_new, rho_a, rho_b,log=TRUE)-dgamma(rho_val, rho_a, rho_b,log=TRUE)
        
        if(runif(1)<exp(lograte)){ ## If accept
            rho <- rho_new
            rho_val <- rho_val_new
        }
    }

    list(rho=rho,rho_val=rho_val)

}

## 
update_bthometies <- function(wins,losses,ties,indices1,indices2,T,K,V,lambda,alpha,theta,at,b){

    ## Sample Z given lambda
    Z1 <- Z2 <- vector("list",T)
    bt <- matrix(0,nrow=K,ncol=T)

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
        
        Z1[[t]] <- sparseMatrix(i=i1, j=j1, x=rgamma(N1,shape=val1, rate=(alpha*lambda[i1, t]+theta*lambda[j1,t])), dims=c(K, K),check=FALSE)
        Z2[[t]] <- sparseMatrix(i=i2, j=j2, x=rgamma(N2,shape=val2, rate=(alpha*lambda[i2, t]+theta*lambda[j2,t])), dims=c(K, K),check=FALSE)
        
        ## check this
        bt[,t] <- b + rowSums(Z1[[t]]) + theta*colSums(Z1[[t]]) + alpha*theta*rowSums(Z2[[t]]) + colSums(Z2[[t]])

    }

    ## Sample lambda given Z,U,theta,alpha
    if(T==1){
        lambda[,1] <- rgamma(K,at[,1],bt[,1])
    }
    else{
        lambda[,1] <- rgamma(K,at[,1]+V[,1],bt[,1]+rho[,1])
        lambda[,T] <- rgamma(K,at[,T]+V[,T-1],bt[,T]+rho[,T-1])
        lambda[,2:(T-1)] <- rgamma(K*(T-2),at[,2:(T-1)]+V[,2:(T-1)]+V[,1:(T-2)],
                                   bt[,2:(T-1)]+rho[,2:(T-1)]+rho[,1:(T-2)])
    }
    
    ## sample alpha given Z,U,lambda
    b_alpha <- 0
    for(t in 1:T){
        b_alpha <- b_alpha + lambda[, t]%*%(rowSums(Z1[[t]]) + theta*rowSums(Z2[[t]]))
    }

    alpha <- rgamma(1,total_wins_and_ties,b_alpha)
    
    
    ## sample theta given Z,U,lambda with MH step
    theta_new <- rnorm(1,theta,.1) # proposal ## use truncated normal to improve eff
    if( theta_new>1){
        temp <- 0
        for( t in 1:T)
            temp <- temp + lambda[, t]%*%(colSums(Z1[[t]]) + alpha*colSums(Z2[[t]]))

        lograte <-  total_ties*log(theta_new^2-1) - theta_new*temp-total_ties*log(theta^2 - 1+.Machine$double.eps) + theta*temp
        if(runif(1)<exp(lograte)) # If accept
            theta <- theta_new
    }

    list(lambda=lambda,alpha=alpha,theta=theta)
}

