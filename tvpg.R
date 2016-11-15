library(BayesLogit)

tvpg <- function(wins, losses, ties, a, bmask.list=NULL,
                          kappa=1, nu=1,omega=Inf,
                          rho, rho_param, Ngibbs, Nburn=1,
                          init.list=NULL,sample.skip=10,drawRho=FALSE){


     ## Number of time steps, T, and number of individuals K
    T <- length(wins)
    K <- nrow(wins[[1]])


    if(is.null(bmask.list)){
        bmask.list[["ALL"]] <- matrix(TRUE,nrow=K,ncol=T,dimnames=list(rownames(wins[[1]]),1:T))
    }
        
    ## Initialization
    if(is.null(init.list)){
        init.list <- init_bthometies(wins,losses,ties,a,bmask.list)
    }

    ## lambda <- init.list$lambda
    ## indices1 <- init.list$indices1
    ## indices2 <- init.list$indices2
    ## theta <- init.list$theta
    ## alpha <- init.list$alpha
    ## at <- init.list$at
    ## bi <- init.list$bi
    ## bmat <- init.list$bmat
    ## b.groups <- matrix(0,nrow=length(bmask.list),ncol=T,dimnames=list(names(bmask.list),1:T))


    ## Store every skip-th value
    Nsamples <- (Ngibbs-Nburn)/sample.skip
    theta.samps <- rep(0,Nsamples)
    alpha.samps <- rep(0,Nsamples)
    beta.samps <- array(0,dim=c(K,T,Nsamples))
    
    lambda.samps[,,1] <- matrix(0,nrow=K,ncol=T)
    for(t in 1:T){
        lambda.samps[,t,1] <- lambda[,t]/sum(lambda[,t])*K 
    }

    theta.samps[1] <- init.list$theta
    alpha.samps[1] <- init.list$alpha
    a.samps[1] <- a
    b.samps[,,1] <- bmat

    ## Run sampler for NGibbs iterations
    for(i in 1:Ngibbs){


        ## TODO: SAMPLING!!
        
    }


}



update_bthometies <- function(wins,losses,ties,indices1,indices2,T,K,V,lambda,alpha,theta,at,bi){

    ## Sample Z given lambda
    Z1 <- Z2 <- vector("list",T)
    bdat <- matrix(0,nrow=K,ncol=T)

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
        
        Z1[[t]] <- sparseMatrix(i=i1, j=j1, x=rpg(N1,val1,beta[i]-beta[j]+log(alpha)+log(theta)), dims=c(K, K),check=FALSE)
        Z2[[t]] <- sparseMatrix(i=i2, j=j2, x=rpg(N2,val2,beta[i]-beta[j]+log(alpha)+log(theta)), dims=c(K, K),check=FALSE)
        
        ## TODO check this and update to include bi's
        bdat[,t] <- rowSums(Z1[[t]]) + theta*colSums(Z1[[t]]) + alpha*theta*rowSums(Z2[[t]]) + colSums(Z2[[t]])

    }

   

    list(lambda=lambda,alpha=alpha,theta=theta)
}

