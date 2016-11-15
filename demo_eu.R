
loadSoccerData <- function() {

    initDir <- getwd()
    setwd("~/Dropbox/tvranking/Data/Soccer/Football-Data.co.uk/")
    years <- list.files()

    tab.full <- c()
    leagues.map <- list()
    for(year in years){
        print(year)
        leagues <- list.files(year,pattern="csv")
        for(league in leagues){
            print(league)
            tab <- read.csv(sprintf("%s/%s",year,league),stringsAsFactors=FALSE)
            tab <- tab[,c("Div","Date","HomeTeam","AwayTeam","FTR")]
            tab.full <- rbind(tab.full,tab)

            tms <- unique(c(tab$HomeTeam,tab$AwayTeam))
            for(tm in tms){
                y <- unlist(strsplit(year,":"))[1]
                leagues.map[[tm]][[y]] <- unlist(strsplit(league,split=".",fixed=TRUE))[1]
            }
        }
    }

    setwd("~/Dropbox/tvranking/Data/Soccer/UEFA/")
    tab.uefa <- c()
    tourneys <- list.files()
    for(tourney in tourneys){
        print(tourney)
        tab <- read.csv(tourney,stringsAsFactors=FALSE)
        tab.uefa <- rbind(tab.uefa,tab)
    }
    names.map <- read.csv("~/Dropbox/tvranking/Data/Soccer/names-map.csv",comment.char="#",stringsAsFactors=FALSE,header=FALSE,allowEscapes=TRUE)

    for(i in 1:nrow(names.map)){
        tab.uefa$Home <- gsub(names.map[i,2],names.map[i,1],tab.uefa$Home)
        tab.uefa$Away <- gsub(names.map[i,2],names.map[i,1],tab.uefa$Away)
    }
    tab.uefa$FTR <- ifelse(tab.uefa$Home.Goals>tab.uefa$Away.Goals,"H","A")
    tab.uefa$FTR[tab.uefa$Home.Goals==tab.uefa$Away.Goals] <- "D"
    tab.uefa$Div <- "UEFA"
    

    Month <- as.numeric(sapply(tab.full$Date,function(x) unlist(strsplit(x, split="/",fixed=TRUE))[2]))
    Year <- as.numeric(sapply(tab.full$Date,function(x) unlist(strsplit(x, split="/",fixed=TRUE))[3]))
    Year <- ifelse(Year<15,Year+2000,Year)
    Year <- ifelse(Year<100,Year+1900,Year)
    Day <- as.numeric(sapply(tab.full$Date,function(x) unlist(strsplit(x, split="/",fixed=TRUE))[1]))

    tab.full <- cbind(tab.full,Day,Month,Year)
    colnames(tab.full)[3] <- "Home"
    colnames(tab.full)[4] <- "Away"

    tab.full <- tab.full[!is.na(tab.full$Year),]

    cols <- c("Home","Away","FTR","Div","Day","Month","Year")
    tab.all <- rbind(tab.full[,cols], tab.uefa[,cols])
    tab.all <- tab.all[order(tab.all$Year,tab.all$Month,tab.all$Day),]

    all.teams1 <- union(unique(tab.full$Home),unique(tab.full$Away))
    all.teams2 <- union(unique(tab.uefa$Home),unique(tab.uefa$Away))
    intersect.teams <- intersect(all.teams1,all.teams2)

    tab.intersect <- tab.all[tab.all$Home%in%intersect.teams&tab.all$Away%in%intersect.teams,]

    leagues.map <- leagues.map[intersect.teams]
    leagues.list <- sort(unique(unlist(leagues.map)))

    setwd(initDir)

    init <- init_wlt(tab.intersect, leagues.list, leagues.map)
    init$leagues.map <- leagues.map
    init$leagues.list <- leagues.list
    
    init
}


## dt in months
init_wlt <- function(tab, leagues.list, leagues.map, dt=1){

    teams <- sort(union(unique(tab$Home),unique(tab$Away)))
    dimnames <- list(teams,teams)
    K <- length(teams)
    nyears <- max(tab$Year)-1994+1
    
    bmask.list <- list()
    for(league in leagues.list){
        bmask.list[[league]] <- matrix(FALSE,nrow=K,ncol=12*nyears/dt,dimnames=list(teams,1:(12*nyears/dt)))
    }
    
    library(Matrix)
    wins <- losses <- ties <- list()
    count <- 1
    for(y in 1994:max(tab$Year)){
        for(m in seq(1,12,by=dt)){
            tab.cur <- tab[tab$Year==y,]
            tab.cur <- tab.cur[tab.cur$Month>=m&tab.cur$Month<(m+dt),]

            w.indices <- which(tab.cur$FTR=="H")
            l.indices <- which(tab.cur$FTR=="A")
            t.indices <- which(tab.cur$FTR=="D")

            ## CHeck for duplicate wins!!
            ## which(duplicatedpaste(tab.cur$Home,tab.cur$Away,tab.cur$FTR)
            
            wins.t <- sparseMatrix(
                i=match(tab.cur$Home[w.indices],teams),
                j=match(tab.cur$Away[w.indices],teams),
                dims=c(K,K),dimnames=dimnames)
            losses.t <- sparseMatrix(
                i=match(tab.cur$Home[l.indices],teams),
                j=match(tab.cur$Away[l.indices],teams),
                dims=c(K,K),dimnames=dimnames)
            ties.t <- sparseMatrix(
                i=match(tab.cur$Home[t.indices],teams),
                j=match(tab.cur$Away[t.indices],teams),
                dims=c(K,K),dimnames=dimnames)

            wins[[count]] <- wins.t
            losses[[count]] <- losses.t
            ties[[count]] <- ties.t
            
            count <- count+1
        }
        
    }
    
    for(tm in teams){
        tm.divs <- rep("",nyears)
        names(tm.divs) <- 1994:max(tab$Year)
        tm.divs[names(leagues.map[[tm]])] <- leagues.map[[tm]]
        year <- 1
        for(league in tm.divs){
            if(league==""){
                if(year==1){
                    tm.divs[year] <- tm.divs[match(TRUE,tm.divs!="")]
                }else{
                    tm.divs[year] <- tm.divs[year-1]
                }
                league <- tm.divs[year]
            }
            bmask.list[[league]][tm,(12/dt*(year-1)+1):(12/dt*year)] <- TRUE
            year <- year+1
        }
    }
    
    return(list(T=nyears*12/dt,wins=wins,losses=losses,ties=ties,K=K,bmask.list=bmask.list,dimnames=dimnames))

}

setwd("~/Dropbox/tvranking/")
init <- loadSoccerData()
wins <- init$wins
losses <- init$losses
ties <- init$ties
bmask.list <- init$bmask.list

if( FALSE ) {
###
    source("/Users/afranks/Dropbox/tvranking/R\ Implementation/tvranking/btemhometies.R")
    source("/Users/afranks/Dropbox/tvranking/R\ Implementation/tvranking/tvbtgibbshome.R")

    rho <- 50*matrix(1,nrow=nrow(wins[[1]]),ncol=(length(wins)-1))
    rho_param <- matrix(c(2000,2),nrow=1)
    a <- 2
    init.list <- init_bthometies(wins, losses, ties, a, bmask.list)

    ## Rprof()
    ## load("~/Dropbox/tvranking/Data/Soccer/eu_constant_b.RData")
    ## init.list$lambda <- apply(results$lambda.samps,c(1,2),mean)
    for(league in rownames(init$bmat) ){
        init.list$bi[bmask.list[[league]]] <- results$b.samps[league,1,90]
        init.list$bmat[league,] <- results$b.samps[league,1,90]
    }

    ## No groups
    nogroups.bmask <- bmask.list
    nogroups.bmask[[1]][] <- TRUE
    for(i in 2:length(nogroups.bmask)) {
        nogroups.bmask[[i]][] <- FALSE
    }
    init.list2 <- init.list
    init.list2$bmat[1, ] <- 1
    init.list2$bmat[2:nrow(init.list2$bmat), ] <- 1
    init.list2$bi[] <- 10
    results <- tvbtgibbshome(wins, losses, ties, a=10, nogroups.bmask,
                             kappa=10, nu=10, omega=Inf, rho, rho_param,
                             Ngibbs=1000, Nburn=10, init.list=init.list2,
                             sample.skip <- 10, drawRho=FALSE)
    save(results, file="~/Dropbox/tvranking/Data/Soccer/eu_nogroups_2015.RData")

    ## NO DATA
    wins2 <- wins; ties2 <- ties; losses2 <- losses;
    for(i in 1:length(wins)) {
        wins2[[i]][] <- 0
        losses2[[i]][] <- 0
        ties2[[i]][] <- 0
    }
    init.list3 <- init.list
    init.list3$bi[] <- 1
    init.list3$bmat[] <- 1
    init.list3$lambda[] <- a/init.list$bi
    results <- tvbtgibbshome(wins2, losses2, ties2, a, bmask.list,
                             kappa=10, nu=10,
                             omega=Inf, rho, rho_param, Ngibbs=1000,
                             Nburn=10, sample.skip=10,
                             init.list=init.list3, drawRho=FALSE)
    

    
    ## Run with non-time varying means
    results <- tvbtgibbshome(wins, losses, ties, a, bmask.list, kappa=10, nu=10,
                             omega=Inf,rho, rho_param, Ngibbs=10000, Nburn=10, sample.skip=10,
                             init.list=init.list, drawRho=FALSE)

    sort(rowMeans(2/results$b.samps[,init$T,]))
    results$lambda.samps[102, , 9]
    save(results,file="~/Dropbox/tvranking/Data/Soccer/eu_constant_b_2015.RData")
    ## Rprof(NULL)

    ## Run with time varying means
    init.list$lambda <- results$lambda.samps[, , 90]
    results <- tvbtgibbshome(wins, losses, ties, a, bmask.list, kappa=10, nu=10, omega=500, rho, rho_param, Ngibbs=10000, Nburn=1000, init.list=init.list, drawRho=FALSE, sample.skip=100)
    save(results,file="~/Dropbox/tvranking/Data/Soccer/eu_varying_b_2015.RData")


    init.list$lambda <- results$lambda.samps[, , dim(results$lambda.samps)[3]]
    init.list$alpha <- results$alpha.samps[length(results$alpha.samps)]
    init.list$theta <- results$theta.samps[length(results$theta.samps)]
    init.list$bmat <- results$b.samps[, , dim(results$b.samps)[3]]

    for(league in rownames(init$bmat) ){
        init.list$bi[bmask.list[[league]]] <- results$b.samps[league, , dim(results$b.samps)[3]]
    }
    
    ## Check mean team skills

    ## first time
    sort(a/results$b.samps[,1,99])
    ## mid time
    sort(a/results$b.samps[,120,99])
    ## last time
    sort(a/results$b.samps[,240,99])

    par(mfrow=c(3, 2))
    for(nm in c("E0", "E3" , "EC")) {
        hist(a/results$b.samps[nm, 264, ], breaks=15, main=sprintf("Mean %s", nm), xlim=c(0, 10))
        hist(sqrt(a/results$b.samps[nm, 264, ]^2), breaks=15, main=sprintf("SD %s", nm), xlim=c(0, 5))
    }
    
    lambda.mean <- apply(results$lambda.samps,c(1,2),mean)
    tms <- sample(1:nrow(lambda.mean),10)
    tms <- c(15,131,135,125,32,132,5,100,99,35,17)

    pdf("tv_plot.pdf")
    plotLambdas(lambda.mean[tms,],rownames(wins[[1]])[tms])
    dev.off()

    german.teams <- which(bmask.list[["D1"]][,240])
    mean(lambda.mean[german.teams,240])
    plotLambdas(lambda.mean[german.teams,],rownames(wins[[1]])[german.teams])

    spanish.teams <- which(bmask.list[["SC0"]][,220])
    mean(lambda.mean[spanish.teams,240])
    plotLambdas(lambda.mean[spanish.teams,],rownames(wins[[1]])[spanish.teams])

    english.teams <- which(bmask.list[["E0"]][,1])
    mean(lambda.mean[english.teams,240])
    plotLambdas(lambda.mean[english.teams,],rownames(wins[[1]])[english.teams])

    plotLambdas(lambda.mean,rownames(wins[[1]]))

    par(mfrow=c(1,1))
    mean.b <- apply(results$b.samps,c(1,2),mean)
    matplot(a/t(mean.b[c("EC", "E0","F1","SC0","SP1","D1","I1","N1"),]),type="l")
    legend("bottomleft",c("EC", "E0","F1","SC0","SP1","D1","I1","N1"),col=1:6,lty=1:5)


    plot(2/mean.b["E0",],type="l",lwd=2,ylim=c(20,70))
    lines(colMeans(lambda.mean*bmask.list[["E0"]]),type="l",lwd=1,col="red")

    lines(mean.b["D1",],type="l",ylim=c(0,10),col="green")
    lines(mean.b["E3",],ylim=c(0,10),col="blue")
    lines(mean.b["SP1",],ylim=c(0,10),col="red")

    leagues <- names(bmask.list)
    leagues <- c("E0","P1","SP1","D1","F1","SC0","I1")
    plot(colMeans(lambda.mean*bmask.list[["E0"]]),type="l",ylim=c(0,0.2),lwd=0)
    count <- 1
    ltys <- rep(1:5,length.out=length(leagues))
    cols=rep(1:6,length.out=length(leagues))
    for(nm in leagues){
        lines(colMeans(lambda.mean*bmask.list[[nm]]),type="l",col=cols[count],lty=ltys[count])
        count <- count+1
    }

    legend("topright",leagues,lty=rep(1:5,length.out=length(leagues)),col=rep(1:6,length.out=length(leagues)),cex=.7,lwd=2)


    tmp <- apply(results$b.samps,c(1,2),mean)
    plotLambdas(tmp,colnames(tmp))


}
