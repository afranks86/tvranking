gen_tvbthometies <- function(K,T,a,b,rho,alpha,theta,lambda=NULL,games=2){
    library(Matrix)
    wins <- losses <- ties <- list()
    for(t in 1:T){
        wins[[t]] <- losses[[t]] <- ties[[t]] <- Matrix(0,nrow=K,ncol=K)
    }

    drawLambda <- FALSE
    if(is.null(lambda)){
        ## Initialize skill
        lambda <- matrix(0,nrow=K,ncol=T)
        lambda[,1] <- rgamma(K,a,b)
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
        
        ## Update lambda's
        if(drawLambda & t<T){
            V <- rpois(K,rho[,t]*lambda[,t])
            lambda[,t+1] <- rgamma(K,a+V,rho[,t]+b)
        }

    }

    list(wins=wins,ties=ties,losses=losses,lambda=lambda)
}

if(FALSE){

    K <- 3
    T <- 100
    a <- 2
    b <- 2
    rho <- matrix(250,nrow=K,ncol=(T-1))
    alpha <- 1
    theta <- 1

    data <- gen_tvbthometies(K,T,a,b,rho,alpha,theta,lambda=Matrix(data=rep(c(0.5,1,2),100),nrow=3),games=20)
    

    ### Advanced Test
    source("tvbtgibbshome.R")
    source("btemhometies.R")
    library(ggplot2)
    
    K <- 10
    T <- 100
    a <- 2
    b <- 2

    rho_param <- matrix(c(800,2,1,.1),nrow=2,byrow=TRUE)
    
    rho <- matrix(rgamma(1,rho_param[1,1],rho_param[1,2]),nrow=K,ncol=(T-1))
    rho[,c(25,50,75)] <- rgamma(1,rho_param[2,1],rho_param[2,2])
    alpha <- 1.5
    theta <- 3 
    
    data <- gen_tvbthometies(K,T,a,b,rho,alpha,theta,lambda=NULL,games=20)

    init.list <- init_bthometies(data$wins,data$losses,data$ties,a)

    rho2 <- matrix(rgamma(1,rho_param[1,1],rho_param[1,2]),nrow=K,ncol=(T-1))
    rho2[,c(25,50,75)] <- rgamma(1,rho_param[2,1],rho_param[2,2])
    results <- tvbtgibbshome(data$wins,data$losses,data$ties,a,rho2,rho_param,Ngibbs=1000,Nburn=100,init.list=init.list,drawRho=TRUE)

    df <- as.data.frame(list(x=rep(1:100,K),lambda=as.vector(1.5*t(data$lambda)),team=rep(1:K,each=T)))
    df$lower.est <- as.vector(t(apply(results$lambda,c(1,2),function(x) quantile(x,0.025))))
    df$upper.est <- as.vector(t(apply(results$lambda,c(1,2),function(x) quantile(x,0.975))))
    df$mean.est <- as.vector(t(apply(results$lambda,c(1,2),mean)))

    pdf("simulation.pdf")
    ggplot(data=df)+geom_line(aes(x=x,y=lambda))+facet_wrap(~team)+geom_ribbon(aes(x=x,ymin=lower.est,ymax=upper.est),alpha=0.3)+geom_line(aes(x=x,y=mean.est),colour="red")
    dev.off()
    
    ## ggplot()+geom_ribbon(aes(x=1:100,ymin=quantile(results$theta,0.025),ymax=quantile(results$theta,0.975)),alpha=0.3)

    ## ggplot()+geom_ribbon(aes(x=1:100,ymin=quantile(results$alpha,0.025),ymax=quantile(results$alpha,0.975)),alpha=0.3)


    ### Debug sampleRho
    rho_mask <- array(0,dim=c(K,T-1,1))
    rho_mask[,,1] <- rho[1]
    sampleRho(data$lambda,rho,rho[1],rho_mask,a,b,rho_param[1],rho_param[2])
    
}
