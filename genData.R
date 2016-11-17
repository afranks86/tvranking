library(Matrix)

gen_tvbthometies <- function(K, T, a, b=NULL, bmask.list, rho,
                             alpha, theta, lambda=NULL,
                             nu=1, kappa=1, omega=Inf, games=2) {

    wins <- losses <- ties <- list()
    for(t in 1:T){
        wins[[t]] <- losses[[t]] <- ties[[t]] <- Matrix(0,nrow=K,ncol=K)
    }

    drawB <- FALSE
    if(is.null(b)){
        ## initialize group means
        b <- matrix(0, dim(bmask.list))
        b[, 1] <- rgamma(nrow(bmask.list), kappa, nu)
        drawB <- TRUE
    }

    drawLambda <- FALSE
    if(is.null(lambda)){
        ## Initialize skill
        lambda <- matrix(0, nrow=K, ncol=T)
        lambda[, 1] <- rgamma(K, a, b[, 1])
        V <- matrix(0, nrow=K, ncol=(T-1))
        drawLambda <- TRUE
    }

    
    for(t in 1:T){
        ## Generate w,t,l

        ## random matchups
        for(i in 1:games){
            home <- sample(1:K,K/2)
            away <- sample(setdiff(1:K,home),K/2)

            home.win.prob <- alpha*lambda[home,t]/(alpha*lambda[home,t]+theta*lambda[away,t])
            away.win.prob <- lambda[away,t]/(alpha*theta*lambda[home,t]+lambda[away,t])

            outcome <- apply(runif(floor(K/2))<cbind(home.win.prob,home.win.prob+away.win.prob),1,sum)
            win.indices <- which(outcome==2)
            loss.indices <- which(outcome==1)
            tie.indices <- which(outcome==0)

            wins[[t]] <- wins[[t]]+sparseMatrix(i=home[win.indices],j=away[win.indices],dims=c(K,K))
            losses[[t]] <- losses[[t]]+sparseMatrix(i=home[loss.indices],j=away[loss.indices],dims=c(K,K))
            ties[[t]] <- ties[[t]]+sparseMatrix(i=home[tie.indices],j=away[tie.indices],dims=c(K,K))

        }

        ## Update b's
        if(drawB & t<T){
            W <- rpois(nrow(b), nu*omega*b[,t])
            b[, t+1] <- rgamma(nrow(b),kappa+W,nu*(1+omega))
        }
        ## Update lambda's
        if(drawLambda & t < T){
            V[, t] <- rpois(K, b[, t] * rho[, t] * lambda[,t])
            lambda[, t+1] <- rgamma(K, a + V[, t], b[, t+1] * (1 + rho[, t]))
        }

    }

    list(wins=wins, ties=ties, losses=losses,
         lambda=lambda, V=V, b=b)
}

if(FALSE){
    ##############################
    ## Single group simulation
    ##############################
    source("tvbtgibbshome.R")
    source("btemhometies.R")
    library(ggplot2)
    
    K <- 10
    T <- 100
    a <- 2

    rho_param <- matrix(c(100, 2, 1, .1), nrow=2, byrow=TRUE)
    
    rho <- matrix(rgamma(1, rho_param[1, 1], rho_param[1,2]), nrow=K, ncol=(T-1))
    ## rho[, c(25, 50, 75)] <- rgamma(1, rho_param[2, 1], rho_param[2, 2])
    alpha <- 1.5
    theta <- 3

    b <- matrix(2, nrow=1, ncol=T)
    bmask.list = list("group"=matrix(TRUE, nrow=K, ncol=T))

    
    
    data <- gen_tvbthometies(K=K, T=T, a=a, b=b, bmask.list=bmask.list,
                             rho=rho, alpha=alpha, theta=theta,
                             lambda=NULL, games=20)
    df <- as.data.frame(list(x=rep(1:100,K),lambda=as.vector(t(data$lambda)),team=rep(1:K,each=T)))
    ggplot(data=df) + geom_line(aes(x=x,y=lambda)) + facet_wrap(~team)
    
    ## verify mean is a/b
    mean(rowMeans(data$lambda))
    
    ## verify var is approx a/b^2
    var(as.numeric(data$lambda))

    ## verify regression to mean: should be close to a/b
    mean(sapply(1:K, function(k) {
    (rho[1, 1] + 1)*mean(data$lambda[k, 2:T] - rho[1, 1]/(1 + rho[1,1])*data$lambda[k, 1:(T-1)]) }))
    
    init.list <- init_bthometies(data$wins, data$losses, data$ties, a, bmask.list)

    results <- tvbtgibbshome(data$wins, data$losses, data$ties,
                             a=a, bmask.list=bmask.list,
                             rho=rho, rho_param=rho_param, Ngibbs=1000,
                             Nburn=100, init.list=init.list,
                             draw=c(b=FALSE, V=TRUE, lambda=TRUE))


    ## rescale data
    nsamps <- dim(results$lambda.samps)[3]
    for(i in 1:nsamps){
        rescaler <- colSums(data$lambda) / colSums(results$lambda.samps[, , i])
        results$lambda.samps[, , i] <- t(t(results$lambda.samps[, , i])*rescaler)
    }
    mean.est <- apply(results$lambda.samps, c(1,2), mean)

    df <- as.data.frame(list(x=rep(1:100,K),lambda=as.vector(t(data$lambda)),team=rep(1:K,each=T)))
    df$lower.est <- as.vector(t(apply(results$lambda,c(1,2),function(x) quantile(x,0.025))))
    df$upper.est <- as.vector(t(apply(results$lambda,c(1,2),function(x) quantile(x,0.975))))
    df$mean.est <- as.vector(t(apply(results$lambda,c(1,2),mean)))

    1 - (mean(df$upper.est < df$lambda) +  mean(df$lower.est > df$lambda))
    
    pdf("single_simulation.pdf")

    ggplot(data=df) + geom_line(aes(x=x,y=lambda)) +
        facet_wrap(~team) +
        geom_ribbon(aes(x=x,ymin=lower.est,ymax=upper.est),alpha=0.3) +
        geom_line(aes(x=x,y=mean.est),colour="red")
    
    dev.off()



    ##############################
    ## Three group mean variation
    ##############################

    K <- 30
    T <- 100
    a <- 2

    rho_param <- matrix(c(100, 2, 1, .1), nrow=2, byrow=TRUE)
    
    rho <- matrix(rgamma(1, rho_param[1, 1], rho_param[1,2]), nrow=K, ncol=(T-1))
    alpha <- 1.5
    theta <- 3

    b <- matrix(2, nrow=3, ncol=T)
    bmask.list = list("group"=matrix(TRUE, nrow=K, ncol=T))
        
    
}
