http://football-data.co.uk/mmz4281/1516/E0.csv

for(year in c("1415","1516")){
    for(league in leagues2) {
        write.csv(read.csv(sprintf("http://football-data.co.uk/mmz4281/%s/%s",year,league)), file=sprintf("~/Dropbox/tvranking/Data/Soccer/Football-Data.co.uk/%s/%s",year,league))
    }
}
