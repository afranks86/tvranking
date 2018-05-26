rm(list=ls())

## Load Soccer Data
source("btemhometies.R")
source("tvbtgibbshome.R")
source("demo_eu.R")
source("holdOut.R")

init <- loadSoccerData()

testTrainSplit <- RemoveRandomGames(init$wins, init$losses, init$ties, init$bmask.list, fracToRemove=0.05, interleague=FALSE)

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
                             a=1, init$bmask.list)

rho <- 200*matrix(1,nrow=nrow(winsTrain[[1]]),ncol=(length(winsTrain)-1))
rho_param <- matrix(c(2000,2),nrow=1)
a <- 1

for(league in rownames(init$bmat) ){
    init.list$bi[bmask.list[[league]]] <- results$b.samps[league,1,90]
    init.list$bmat[league,] <- results$b.samps[league,1,90]
}


## Init List
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
init.list2$bmat <- matrix(1, nrow=1, ncol=init.list2$T,
                          dimnames=list(c("group 1"), 1:init.list2$T))
init.list2$bi[] <- 1
init.list2$lambda[] <- 1
results <- tvbtgibbshome(winsTrain, lossesTrain, tiesTrain,
                         a=2, bmask.list=nogroups.bmask, kappa=1, nu=1,
                         omega=Inf, rho, rho_param, Ngibbs=1000,
                         Nburn=100, init.list=init.list2, draw=c(b=FALSE))

ll.nogroups <- HeldOutLikelihood(init$K, winsTest, lossesTest, tiesTest, results$lambda.samps, results$alpha.samps, results$theta.samps, meanOnly=TRUE)

lamsame <- results$lambda.samps
lamsame[] <- 1/2
asame <- results$alpha.samps
## asame[] <- 1

tsame <- results$theta.samps
## tsame[] <- 5

ll.even <- HeldOutLikelihood(init$K, winsTest, lossesTest, tiesTest, lamsame, asame, tsame, meanOnly=TRUE)

## Run with groups, non-time varying means
resultsGroups <- tvbtgibbshome(winsTrain, lossesTrain, tiesTrain, a, init$bmask.list,
                         kappa=1, nu=1, omega=Inf, rho, rho_param, Ngibbs=1000, Nburn=100, init.list=init.list, drawRho=FALSE)

ll.groups <- HeldOutLikelihood(init$K, winsTest, lossesTest, tiesTest, resultsGroups$lambda.samps, resultsGroups$alpha.samps, resultsGroups$theta.samps, meanOnly=TRUE)

## Run with time varying means
results <- tvbtgibbshome(winsTrain, lossesTrain, tiesTrain, a, init$bmask.list, kappa=1, nu=1, omega=500, rho, rho_param, Ngibbs=100, Nburn=10, init.list=init.list, drawRho=FALSE, sample.skip=10)
n
ll.tvmeans <- HeldOutLikelihood(init$K, winsTest, lossesTest, tiesTest, results$lambda.samps, results$alpha.samps, results$theta.samps)




## first time
sort(results$b.samps[,1,9])
## mid time
sort(results$b.samps[,132,9])
## last time
sort(results$b.samps[,264,9])

lambda.mean <- apply(results$lambda.samps,c(1,2),mean)
lambda.mean <- results$lambda.samps[, , 90]

tms <- c(15,131,135,125,32,132,5,100,99,35,17)

tms <- sample(1:nrow(lambda.mean),10)
plotLambdas(lambda.mean[tms,],rownames(wins[[1]])[tms])
