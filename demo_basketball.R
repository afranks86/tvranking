### Rows are home teams
files <- list.files("~/course/NBA_pbp/")

wins <- losses <- ties <- list()

teams <- c("Atlanta","Boston","Brooklyn","Charlotte","Chicago","Cleveland",
           "Dallas","Denver","Detroit","Golden State","Houston","Indiana",
           "LA Clippers","LA Lakers","Memphis","Miami","Milwaukee","Minnesota",
           "New Orleans","New York","Oklahoma City","Orlando","Philadelphia",
           "Phoenix","Portland","Sacramento","San Antonio","Toronto","Utah",
           "Washington")
K <- 30
T <- 48

## Init list
for(t in 1:48){
    wins[[t]] <- losses[[t]] <- ties[[t]] <- Matrix(0,nrow=K,ncol=K,dimnames=list(rwnms=teams,colnms=teams))
}

count <- 1
for(file in files){
    print(count)

    game.table <- read.csv(sprintf("~/course/NBA_pbp/%s",file),check.names=FALSE)
    
    home <- colnames(game.table)[3]
    away <- colnames(game.table)[2]

    indices <- 49-findInterval(game.table$Time,seq(0,2880,by=60))
    last.idx.in.interval <- length(indices)+1-match(1:48,rev(indices))
    last.score <- c(0,0)
    for(t in 1:48){
        idx <- last.idx.in.interval[t]
        if(is.na(idx))
            ties[[t]][home,away] <- ties[[t]][home,away]+1
        else{
            score.change <- game.table[idx,2:3]-last.score
            sgn <- sign(score.change[2]-score.change[1])
            if(sgn>0)
                wins[[t]][home,away] <- wins[[t]][home,away]+1
            else if(sgn<0)
                losses[[t]][home,away] <- losses[[t]][home,away]+1
            else
                ties[[t]][home,away] <- ties[[t]][home,away]+1

            last.score <- game.table[idx,2:3]
        }
    }
    count <- count+1
}

## TODO get wins losses 

source("/Users/afranks/Dropbox/tvranking/R\ Implementation/tvranking/btemhometies.R")
source("/Users/afranks/Dropbox/tvranking/R\ Implementation/tvranking/tvbtgibbshome.R")

rho <- 2000*matrix(1,nrow=nrow(wins[[1]]),ncol=(length(wins)-1))
rho[,c(12,24,36)] <- 90
rho_param <- matrix(c(200,2,10,.1),nrow=2,byrow=TRUE)
a <- 2
init.list <- init_bthometies(wins,losses,ties,a)

#Rprof()
results <- tvbtgibbshome(wins,losses,ties,a,rho,rho_param,Ngibbs=1000,Nburn=100,init.list=init.list,drawRho=FALSE)
#Rprof(NULL)

lambda.mean <- apply(results$lambda.samps,c(1,2),mean)
idx <- match(c("Brooklyn","Chicago","Miami","Indiana","Philadelphia"),rownames(wins[[1]]))
idx <- sample(1:length(rownames(wins[[1]])),10)
plotLambdas(lambda.mean[idx,],rownames(wins[[1]])[idx])

q1 <- rowMeans(lambda.mean[,1:12])
q2 <- rowMeans(lambda.mean[,13:24])
q3 <- rowMeans(lambda.mean[,25:36])
q4 <- rowMeans(lambda.mean[,37:48])

names(q1) <- names(q2) <- names(q3) <- names(q4) <- rownames(wins[[1]])

tmp <- t(t(lambda.mean)/colSums(lambda.mean))
tmp2 <- apply(tmp,1,sd)
names(tmp2) <- names(q1)
sort(tmp2)

