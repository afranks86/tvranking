library(plyr)

tab <- read.csv("../../../Data/Golf/revent.TXT",sep=";",header=TRUE,stringsAsFactors=FALSE,quote="")
tab <- tab[,1:27]
allPlayers <- sort(unique(tab$Player.Name))
allTournaments <- unique(tab$Event.Name)

T <- length(allTournaments)
K <- length(allPlayers)

ranks <- matrix(nrow=K,ncol=T,dimnames=list(allPlayers,allTournaments))

for(tourney in allTournaments) {

    indices <- which(tab$Event.Name==tourney)
    players <- tab$Player.Name[indices]
    tourney.ranks <- tab$Finish.Position.numeric[indices]
    ranks[players,tourney] <- tourney.ranks
    
}

