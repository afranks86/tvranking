init_bthometies <- function(wins, losses, ties, a, bmask.list){

    K <- nrow(wins[[1]])
    T <- length(wins)
    indices1 <- indices2 <- vector("list",T)
    wltStat <- matrix(0,nrow=K,ncol=T)
    bi <- matrix(0,nrow=K,ncol=T,dimnames=list(rownames(wins[[1]]),1:T))
    bmat <- matrix(0,nrow=length(bmask.list),ncol=T,dimnames=list(names(bmask.list),1:T))
    for(t in 1:T){

        wt <- wins[[t]] + ties[[t]]
        lt <- losses[[t]]+ties[[t]]
        
        indices1[[t]] <- summary(wt)
        indices2[[t]] <- summary(lt)

        N1 <- nrow(indices1[[t]])
        i1 <- indices1[[t]]$i
        j1 <- indices1[[t]]$j
        val1 <- indices1[[t]]$x

        N2 <- nrow(indices2[[t]])
        i2 <- indices2[[t]]$i
        j2 <- indices2[[t]]$j
        val2 <- indices2[[t]]$x

        wltStat[, t] <- rowSums(wt+t(lt))
    }
    
    ## Initial values for lambda, theta and alpha
    ## Init lambda with a sliding window EM
    prec <- 1e-5
    lag <- 30
    fit <- swbtemhometies(wins, losses, ties, a, prec, lag)

    
    for(group in names(bmask.list)){
        bmask <- bmask.list[[group]]
        bval <- a / mean(fit$lambda[bmask])
        bmat[group,] <- bval
        bi[bmask] <- bval
    }

    lambda <- fit$lambda*a
    
    list(K=K,T=T, indices1=indices1, indices2=indices2,
         lambda=fit$lambda, alpha=fit$alpha, theta=fit$theta, wltStat=wltStat,
         bi=bi, bmat=bmat)
    
}

swbtemhometies <- function (wins, losses, ties, a, prec=1e-8, lag){

    ## Sliding windows EM for time varying ranking
    ## Run sliding-window EM
    K = nrow(wins[[1]])
    T = length(wins)

    lambda_em <- matrix(0,nrow=K, ncol=T)
    theta_em <- rep(0,T)
    alpha_em <- rep(0,T)
    for(t in 1:T){
        print(t)
        w.EM <- Matrix(0,nrow=K, ncol=K)
        t.EM <- Matrix(0,nrow=K, ncol=K)
        l.EM <- Matrix(0,nrow=K, ncol=K)
        for(m in max(1, t-lag):min(T, t+lag)){
            w.EM <- w.EM + wins[[m]]
            l.EM <- l.EM + losses[[m]]
            t.EM <- t.EM + ties[[m]]
        }
        fit <- btemhometies(w.EM, l.EM, t.EM, a, prec)
        lambda_em[, t] <- fit$lambda
        theta_em[t] <- fit$theta
        alpha_em[t] <- fit$alpha
    }

    list(lambda=lambda_em,theta=mean(theta_em),alpha=mean(alpha_em))
}


btemhometies <- function(wins, losses, ties, a, prec=1e-8){
###########################################################################
    ## function [lambda.samps, theta.samps, alpha.samps, a.samps, stats] = btemhometies(w, l, t,
    ## a, prec)
    ## EM algorithm for Bradley-Terry models with home advantage and ties
##########

    ## INPUTS
    ##
    ## w(i,j) is the number of times i beats j when i is at home
    ## l(i,j) is the number of times i looses to j when i is at home
    ## t(i, j) is the number of times i ties j when i is at home
    ## a is the shape parameter for the gamma prior (default a = 1)
    ## prec is the precision of the EM (default = 1e-8)

##########
    ## OUTPUTS
    ## lambda.samps gives the values of the normalized skills at each iteration
    ## theta.samps gives the values of the ties parameter at each iteration
    ## alpha.samps gives the values of the home advantage parameter at each iteration
    ## a.samps gives the values of the shape parameter at each iteration
    ## stats is a structure with some summary statistics on the parameters
###########################################################################
    ## Reference:
    ## F. Caron and A. Doucet. Efficient Bayesian inference for generalized
    ## Bradley-Terry models. To appear in Journal of Computational and Graphical
    ## Statistics, 2011.
    ##
    ## February 2011
    ## Author: F. Caron (INRIA Bordeaux Sud-Ouest)
    ## Francois.Caron@inria.fr
    ## http://www.math.u-bordeaux1.fr/~fcaron/
###########################################################################

    K <- nrow(wins)
    iter_max <- 5000
    b <- K*a - 1

    lambda <- rep(1,K)

    lambda.samps <- matrix(0,iter_max, K)
    theta.samps <- rep(0,iter_max)
    alpha.samps <- rep(0,iter_max)
    ell <- rep(0,iter_max)

    s1 = wins + ties
    s2 = losses + ties
    ak = a - 1 + rowSums(s1 + t(s2)) 

    indices1 <- summary(s1)
    i1 <- indices1$i
    j1 <- indices1$j
    val1 <- indices1$x

    indices2 <- summary(s2)
    i2 <- indices2$i
    j2 <- indices2$j
    val2 <- indices2$x

    T = sum(sum(ties))
    H = sum(val1) # Number of ties and wins at home

    change <- rep(Inf,length(lambda))

    theta <- 1.5 # Tie parameter
    alpha <- 1; # Home advantage parameter


    iteration <- 1
    lambda.samps[iteration, ] <- lambda/sum(lambda)
    theta.samps[iteration] <- theta
    alpha.samps[iteration] <- alpha

    ## Compute log-posterior to monitor convergence
    temp1 <- sparseMatrix(i=i1, j=j1, x=val1*(log(lambda[i1])-log(alpha*lambda[i1]+theta*lambda[j1])), dims=c(K, K))
    temp2 <- sparseMatrix(i=i2, j=j2, x=val2*(log(lambda[j2])-log(alpha*theta*lambda[i2]+lambda[j2])), dims=c(K, K))

    ell[iteration] <- sum(sum(val1))*log(alpha) + T*log(theta^2-1)+ sum(sum(temp1)) + sum(sum(temp2)) + (a-1)*sum(log(lambda)) - b*sum(lambda)

    while( sqrt(sum(change^2)) > prec  && iteration<iter_max){

        iteration <- iteration + 1
        
        ## E step
        Z <- sparseMatrix(i=i1, j=j1, x=val1/(alpha*lambda[i1]+theta*lambda[j1]), dims=c(K, K))
        U <- sparseMatrix(i=i2, j=j2, x=val2/(alpha*theta*lambda[i2]+lambda[j2]), dims=c(K, K))
        
        ## Maximize w.r.t. lambda given (theta, alpha)
        bk = b + rowSums(Z) + theta*colSums(Z) + alpha*theta*rowSums(U) + colSums(U)
        lambda_new <- ak/bk
        
        change <- lambda_new/sum(lambda_new) - lambda/sum(lambda)
        lambda <- lambda_new

        ## Maximize w.r.t. alpha given (theta, lambda)
        ## Flat prior on alpha
        alpha <- H/(lambda %*% (rowSums(Z) + theta*rowSums(U)))

        ## Maximize w.r.t. theta given (alpha, lambda)
        ## Flat prior on theta 
        q <- lambda%*%(colSums(Z) + alpha*rowSums(U))
        if(T==0)
            theta <- 1
        else
            theta <- T/q * (1+sqrt(1+q^2/T^2))

        ## Store outputs  
        lambda.samps[iteration,] <- lambda/sum(lambda)
        theta.samps[iteration] <- theta
        alpha.samps[iteration] <- alpha

        ## Compute log-posterior to monitor convergence
        temp1 <- sparseMatrix(i=i1, j=j1, x=val1*(log(lambda[i1])-log(alpha*lambda[i1]+theta*lambda[j1])), dims=c(K, K))
        temp2 <- sparseMatrix(i=i2, j=j2, x=val2*(log(lambda[j2])-log(alpha*theta*lambda[i2]+lambda[j2])), dims=c(K, K))
        ell[iteration] <- sum(sum(val1))*log(alpha) + T*log(theta^2-1)
        + sum(sum(temp1)) + sum(sum(temp2)) + (a-1)*sum(log(lambda)) - b*sum(lambda)
    }

    lambda.samps <- lambda.samps[1:iteration,]
    theta.samps <- theta.samps[1:iteration]
    alpha.samps <- alpha.samps[1:iteration]
    ell <- ell[1:iteration]

    lambda <- lambda/sum(lambda)*K


    list(lambda=lambda,theta=theta,alpha=alpha)

}
