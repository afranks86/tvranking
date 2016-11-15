rm(list=ls())
a <- 1
## Load Soccer Data
source("/Users/afranks/Dropbox/tvranking/R\ Implementation/tvranking/btemhometies.R")
source("/Users/afranks/Dropbox/tvranking/R\ Implementation/tvranking/tvbtgibbshome.R")
source("/Users/afranks/Dropbox/tvranking/R\ Implementation/tvranking/demo_eu.R")

source("/Users/afranks/Dropbox/tvranking/R\ Implementation/tvranking/holdOut.R")

init <- loadSoccerData()

testTrainSplit <- RemoveRandomGames(init$wins, init$losses, init$ties, init$bmask.list, fracToRemove=0.05, interleague=TRUE)

winsTrain <- testTrainSplit$winsTrain
lossesTrain <- testTrainSplit$lossesTrain
tiesTrain <- testTrainSplit$tiesTrain

winsTest <- testTrainSplit$winsTest
lossesTest <- testTrainSplit$lossesTest
tiesTest <- testTrainSplit$tiesTest

## non-tv held-out likelihood
## tv, no groups
## tv, with groups
## tv, with tv groups

init.list <- init_bthometies(winsTrain, lossesTrain, tiesTrain,
                             a, init$bmask.list)

rho <- 50*matrix(1,nrow=nrow(winsTrain[[1]]),ncol=(length(winsTrain)-1))
rho_param <- matrix(c(2000,2),nrow=1)
a <- 2

for(league in rownames(init$bmat) ){
    init.list$bi[bmask.list[[league]]] <- results$b.samps[league,1,90]
    init.list$bmat[league,] <- results$b.samps[league,1,90]
}


## Init List
init.list$lambda
tmp <- array(data=0,dim=c(init$K,init$T,90))
for(i in 1:90){
    tmp[, , i] <-  init.list$lambda
}
alpha.samps = rep(init.list$alpha,90)
theta.samps = rep(init.list$theta,90)
ll.em <- HeldOutLikelihood(init$K, winsTest, lossesTest, tiesTest, tmp, alpha.samps, theta.samps)

## No groups
nogroups.bmask <- init$bmask.list
nogroups.bmask[[1]][] <- 1
for(i in 2:length(nogroups.bmask)) {
    nogroups.bmask[[i]][] <- 0
}
init.list2 <- init.list
init.list2$bmat[1,] <- 1
init.list2$bmat[2:nrow(init.list2$bmat),] <- 0
results <- tvbtgibbshome(winsTrain,lossesTrain,tiesTrain, a=5, nogroups.bmask, kappa=1, nu=1, omega=Inf, rho, rho_param, Ngibbs=100, Nburn=10, init.list=init.list2,drawRho=FALSE) 

ll.nogroups <- HeldOutLikelihood(init$K, winsTest, lossesTest, tiesTest, results$lambda.samps, results$alpha.samps, results$theta.samps, meanOnly=TRUE)

## Run with groups, non-time varying means
resultsGroups <- tvbtgibbshome(winsTrain, lossesTrain, tiesTrain, a, init$bmask.list,
                         kappa=1, nu=1, omega=Inf, rho, rho_param, Ngibbs=1000, Nburn=100, init.list=init.list, drawRho=FALSE)

ll.groups <- HeldOutLikelihood(init$K, winsTest, lossesTest, tiesTest, resultsGroups$lambda.samps, resultsGroups$alpha.samps, resultsGroups$theta.samps, meanOnly=TRUE)

## Run with time varying means
results <- tvbtgibbshome(winsTrain, lossesTrain, tiesTrain, a, init$bmask.list, kappa=1, nu=1, omega=500, rho, rho_param, Ngibbs=100, Nburn=10, init.list=init.list, drawRho=FALSE, sample.skip=10)

ll.tvmeans <- HeldOutLikelihood(init$K, winsTest, lossesTest, tiesTest, results$lambda.samps, results$alpha.samps, results$theta.samps)




## first time
sort(results$b.samps[,1,9])
## mid time
sort(results$b.samps[,132,9])
## last time
sort(results$b.samps[,264,9])

lambda.mean <- apply(results$lambda.samps,c(1,2),mean)
tms <- sample(1:nrow(lambda.mean),10)
tms <- c(15,131,135,125,32,132,5,100,99,35,17)

plotLambdas(lambda.mean[tms,],rownames(wins[[1]])[tms])
