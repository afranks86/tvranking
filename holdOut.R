## non-tv held-out likelihood
## tv, no groups
## tv, with groups
## tv, with tv groups

RemoveRandomGames <- function(wins, losses, ties, bmask.list, fracToRemove=0.05, interleague=FALSE) {

    winsTrain <- wins
    lossesTrain <- losses
    tiesTrain <- ties
    
    winsTest <- wins
    lossesTest <- losses
    tiesTest <- ties
    
    for(t in 1:length(wins)) {

        wt <- wins[[t]]

        if( interleague ) {
            for(lg in names(bmask.list)) {
                intraleagueIndices <- which(bmask.list[[lg]][,t])
                if( length(intraleagueIndices) > 0 )
                    wt[intraleagueIndices, intraleagueIndices] <- FALSE
            }
        }
        
        wt.indices <- which(wt)
        samps <- wt.indices[runif(length(wt.indices)) < fracToRemove]
        winsTrain[[t]][samps] <- FALSE

        winsTest[[t]][] <- FALSE
        winsTest[[t]][samps] <- TRUE
    }
    for(t in 1:length(losses)) {

        lt <- losses[[t]]

        if( interleague ) {
            for(lg in names(bmask.list)) {
                intraleagueIndices <- which(bmask.list[[lg]][,t])
                if( length(intraleagueIndices) > 0 )
                    lt[intraleagueIndices, intraleagueIndices] <- FALSE
            }
        }
        
        lt.indices <- which(losses[[t]])
        samps <- lt.indices[runif(length(lt.indices)) < fracToRemove]
        lossesTrain[[t]][samps] <- FALSE
        
        lossesTest[[t]][] <- FALSE
        lossesTest[[t]][samps] <- TRUE
    }
    for(t in 1:length(ties)) {
        tt <- ties[[t]]

        if( interleague ) {
            for(lg in names(bmask.list)) {
                intraleagueIndices <- which(bmask.list[[lg]][,t])
                if( length(intraleagueIndices) > 0 )
                    tt[intraleagueIndices, intraleagueIndices] <- FALSE

            }
        }
        
        tt.indices <- which(ties[[t]])
        samps <- tt.indices[runif(length(tt.indices)) < fracToRemove]
        tiesTrain[[t]][samps] <- FALSE

        tiesTest[[t]][] <- FALSE
        tiesTest[[t]][samps] <- TRUE
    }
    

    list(winsTrain=winsTrain, lossesTrain=lossesTrain, tiesTrain=tiesTrain,
         winsTest=winsTest, lossesTest=lossesTest, tiesTest=tiesTest)
}

HeldOutLikelihood <- function(K, winsTest, lossesTest, tiesTest, Lambda, alphas, thetas, meanOnly=FALSE ) {

    Nsamps <- dim(Lambda)[3]
    ll <- 0
    
    for(t in 1:length(winsTest)) {

        ## Home team are rows

        ## Wins
        indices <- which(winsTest[[t]])
        homeTeams <- (indices-1) %% K + 1
        awayTeams <- floor( (indices - 1) / K)+1
        
        lam1 <- matrix(Lambda[homeTeams, t, ], ncol=Nsamps)
        lam2 <- matrix(Lambda[awayTeams, t, ], ncol=Nsamps)

        if(meanOnly) {

            lam1 <- rowMeans(lam1)
            lam2 <- rowMeans(lam2)
            alpha <- mean(alphas)
            theta <- mean(thetas)
            
            homeWinProb <-  alpha*lam1/(alpha*lam1+theta*lam2)
            
        } else {

            alphaLam1 <- sweep(lam1, 2, alphas, '*')
            thetaLam2 <- sweep(lam2, 2, thetas, '*')
        
            homeWinProb <- rowMeans((alphaLam1)/(alphaLam1 + thetaLam2))
        }
        
        ## Losses
        indices <- which(lossesTest[[t]])
        homeTeams <- (indices-1) %% K + 1
        awayTeams <- floor( (indices - 1) / K)+1

        lam1 <- matrix(Lambda[homeTeams, t, ], ncol=Nsamps)
        lam2 <- matrix(Lambda[awayTeams, t, ], ncol=Nsamps)

        if(meanOnly) {
            
            lam1 <- rowMeans(lam1)
            lam2 <- rowMeans(lam2)
            
            alphaThetaLam1 <- lam1*mean(alphas)*mean(thetas)
            
            homeLossProb <- lam2/(alphaThetaLam1 + lam2)
            
        } else {

            alphaThetaLam1 <- sweep(lam1, 2, alphas*thetas, '*')
            
            homeLossProb <- rowMeans(lam2/(alphaThetaLam1 + lam2))
        }
        ## Ties

        indices <- which(tiesTest[[t]])
        homeTeams <- (indices-1) %% K + 1
        awayTeams <- floor( (indices - 1) / K)+1
                
        lam1 <- matrix(Lambda[homeTeams, t, ], ncol=Nsamps)
        lam2 <- matrix(Lambda[awayTeams, t, ], ncol=Nsamps)

        if(meanOnly) {

            lam1 <- rowMeans(lam1)
            lam2 <- rowMeans(lam2)
            
            alphaLam1 <- lam1*mean(alphas)
            thetaLam2 <- lam2*mean(thetas)
            alphaThetaLam1 <- lam1*mean(thetas)*mean(thetas)

            homeWinProbT <- (alphaLam1)/(alphaLam1 + thetaLam2)
            homeLossProbT <- lam2/(alphaThetaLam1 + lam2)

            tieProb <- 1 - homeWinProbT-homeLossProbT
        } else {
            alphaLam1 <- sweep(lam1, 2, alphas, '*')
            thetaLam2 <- sweep(lam2, 2, thetas, '*')
            alphaThetaLam1 <- sweep(lam1, 2, alphas*thetas, '*')

            homeWinProbT <- rowMeans((alphaLam1)/(alphaLam1 + thetaLam2))
            homeLossProbT <- rowMeans(lam2/(alphaThetaLam1 + lam2))

            tieProb <- 1 - homeWinProbT-homeLossProbT
        }

        ll <- ll + sum(log(homeWinProb)) +  sum(log(homeLossProb)) + sum(log(tieProb))
    }
    
    ll
}

TeamPredictions <- function(team="Manchester United", K, wins, losses, ties, Lambda, alphas, thetas, meanOnly=FALSE ) {

    Nsamps <- dim(Lambda)[3]
    ll <- 0

    teamMask <- matrix(0, K, K)
    idx <- which(rownames(wins[[t]])==team)
    teamMask[, idx] <- 1
    teamMask[idx, ] <- 1
    
    for(t in 1:length(wins)) {

        ## Home team are rows

        ## Wins
        indices <- which(wins[[t]]*teamMask==1)
        homeTeams <- (indices-1) %% K + 1
        awayTeams <- floor( (indices - 1) / K)+1
        
        lam1 <- matrix(Lambda[homeTeams, t, ], ncol=Nsamps)
        lam2 <- matrix(Lambda[awayTeams, t, ], ncol=Nsamps)

        if(meanOnly) {

            lam1 <- rowMeans(lam1)
            lam2 <- rowMeans(lam2)
            alpha <- mean(alphas)
            theta <- mean(thetas)
            
            homeWinProb <-  alpha*lam1/(alpha*lam1+theta*lam2)
            
        } else {

            alphaLam1 <- sweep(lam1, 2, alphas, '*')
            thetaLam2 <- sweep(lam2, 2, thetas, '*')
        
            homeWinProb <- rowMeans((alphaLam1)/(alphaLam1 + thetaLam2))
        }
        
        ## Losses
        indices <- which(losses[[t]]*teamMask==1)
        homeTeams <- (indices-1) %% K + 1
        awayTeams <- floor( (indices - 1) / K)+1

        lam1 <- matrix(Lambda[homeTeams, t, ], ncol=Nsamps)
        lam2 <- matrix(Lambda[awayTeams, t, ], ncol=Nsamps)

        if(meanOnly) {
            
            lam1 <- rowMeans(lam1)
            lam2 <- rowMeans(lam2)
            
            alphaThetaLam1 <- lam1*mean(alphas)*mean(thetas)
            
            homeLossProb <- lam2/(alphaThetaLam1 + lam2)
            
        } else {

            alphaThetaLam1 <- sweep(lam1, 2, alphas*thetas, '*')
            
            homeLossProb <- rowMeans(lam2/(alphaThetaLam1 + lam2))
        }
        ## Ties

        indices <- which(ties[[t]]*teamMask==1)
        homeTeams <- (indices-1) %% K + 1
        awayTeams <- floor( (indices - 1) / K)+1
                
        lam1 <- matrix(Lambda[homeTeams, t, ], ncol=Nsamps)
        lam2 <- matrix(Lambda[awayTeams, t, ], ncol=Nsamps)

        if(meanOnly) {

            lam1 <- rowMeans(lam1)
            lam2 <- rowMeans(lam2)
            
            alphaLam1 <- lam1*mean(alphas)
            thetaLam2 <- lam2*mean(thetas)
            alphaThetaLam1 <- lam1*mean(thetas)*mean(thetas)

            homeWinProbT <- (alphaLam1)/(alphaLam1 + thetaLam2)
            homeLossProbT <- lam2/(alphaThetaLam1 + lam2)

            tieProb <- 1 - homeWinProbT-homeLossProbT
        } else {
            alphaLam1 <- sweep(lam1, 2, alphas, '*')
            thetaLam2 <- sweep(lam2, 2, thetas, '*')
            alphaThetaLam1 <- sweep(lam1, 2, alphas*thetas, '*')

            homeWinProbT <- rowMeans((alphaLam1)/(alphaLam1 + thetaLam2))
            homeLossProbT <- rowMeans(lam2/(alphaThetaLam1 + lam2))

            tieProb <- 1 - homeWinProbT-homeLossProbT
        }
        
        if(length(homeWinProb) > 0)
            print(sprintf("Wins: %s", paste(homeWinProb, collapse=",")))
        if(length(tieProb) > 0)
            print(sprintf("Ties: %s", paste(tieProb, collapse=",")))
        if(length(homeLossProb) > 0)
            print(sprintf("Losses: %s", paste(homeLossProb, collapse=",")))
        if(!all(c(length(homeWinProb)==0, length(homeLossProb)==0, length(tieProb)==0))) {
            print(t)
            print("-----------")
        }
    }
    

}
