setwd("~/Dropbox/tvranking/Data/Soccer/Football-Data.co.uk/")
years <- list.files()

tab.full <- c()
leagues.list <- list()
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
            leagues.list[[tm]][[year]] <- unlist(strsplit(league,split=".",fixed=TRUE))[1]
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

# -> Sunderland
# -> Reading
# -> Swansea
# Liverpool

cols <- c("Home","Away","FTR","Div","Day","Month","Year")
tab.all <- rbind(tab.full[,cols],
      tab.uefa[,cols])
tab.all <- tab.all[order(tab.all$Year,tab.all$Month,tab.all$Day),]

all.teams1 <- union(unique(tab.full$Home),unique(tab.full$Away))
all.teams2 <- union(unique(tab.uefa$Home),unique(tab.uefa$Away))
intersect.teams <- intersect(all.teams1,all.teams2)

tab.intersect <- tab.all[tab.all$Home%in%intersect.teams&tab.all$Away%in%intersect.teams,]

leagues.list <- leagues.list[intersect.teams]


## dt in months
init_wlt <- function(tab,dt=1){

    teams <- sort(union(unique(tab$Home),unique(tab$Away)))
    dimnames <- list(teams,teams)
    K <- length(teams)
    
    library(Matrix)
    wins <- losses <- ties <- list()
    count <- 1
    for(y in 1994:2013){
        for(m in seq(1,12,by=dt)){
            tab.cur <- tab[tab$Year==y,]
            tab.cur <- tab.cur[tab.cur$Month>=m&tab.cur$Month<(m+dt),]

            w.indices <- which(tab.cur$FTR=="H")
            l.indices <- which(tab.cur$FTR=="A")
            t.indices <- which(tab.cur$FTR=="D")

            # CHeck for duplicate wins!!
            #which(duplicatedpaste(tab.cur$Home,tab.cur$Away,tab.cur$FTR)
            
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

    return(list(wins=wins,losses=losses,ties=ties,K=K,dimnames=dimnames))
}

init <- init_wlt(tab.intersect)
wins <- init$wins
losses <- init$losses
ties <- init$ties

###
source("/Users/afranks/Dropbox/tvranking/R\ Implementation/btemhometies.R")
source("/Users/afranks/Dropbox/tvranking/R\ Implementation/tvbtgibbshome.R")

rho <- 300*matrix(1,nrow=nrow(wins[[1]]),ncol=(length(wins)-1))
rho_param <- matrix(c(2000,2),nrow=1)
a <- 2
init.list <- init_bthometies(wins,losses,ties,a)

#Rprof()
results <- tvbtgibbshome(wins,losses,ties,a,rho,rho_param,Ngibbs=1000,Nburn=100,init.list=init.list,drawRho=TRUE)
#Rprof(NULL)

lambda.mean <- apply(results$lambda.samps,c(1,2),mean)
tms <- sample(1:nrow(lambda.mean),10)
tms <- c(tms,32,132,5,100)
plotLambdas(lambda.mean[tms,],rownames(wins[[1]])[tms])

plotLambdas(lambda.mean,rownames(wins[[1]]))

plotLambdas <- function(lambda.mean,teams){

    matplot(t(lambda.mean),type="l",xlim=c(0,ncol(lambda.mean)+20),lwd=2)
    legend("topright",teams,lty=rep(1:5,length.out=length(teams)),col=rep(1:6,length.out=length(teams)),cex=.7,lwd=2)
}
