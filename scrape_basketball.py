from bs4 import BeautifulSoup
import requests
import csv

## First get all game ID's
pbpLinks = []
for month in range(1,12+1):
    for day in range(1,31+1):
        if month in range(1,7):
            year="2014"
        else:
            year="2013"
        print("Day: %i, Month: %i, Year: %s" % (day,month,year))
        req = requests.get("http://www.basketball-reference.com/boxscores/index.cgi?month=%i&day=%i&year=%s" % (month,day,year))
        data = req.text
        soup = BeautifulSoup(data)
        links = soup.findAll('a',{"href" : re.compile(".+pbp.+")})
        for link in links:
            pbpLinks.append(link.get("href"))

## For each game get the pbp score
for game in pbpLinks:
    req  = requests.get("http://www.basketball-reference.com%s" % game)
    data = req.text
    pbpSoup = BeautifulSoup(data)
    tableHeader = pbpSoup.find('th',text="Time").findParent().findChildren('th')
    away = tableHeader[1].string
    home = tableHeader[3].string

    
    game_id = re.search('/boxscores/pbp/(.+)\.html',game).group(1)
    file = open("/Users/afranks/Dropbox/course/NBA_pbp/%s-%s-%s.csv" % (away,home,game_id), 'w')
    writer = csv.writer(file)
    res = pbpSoup.findAll('td',{"class" : "align_right background_lime"})
    writer.writerow(["Time",away,home])

    time_last = 720.1
    quarter = 1
    for r in res:
        parent = r.findParent()
        children = parent.findChildren('td')

        ## Find time in seconds, adjust for quarter
        time = children[0].string
        times = re.split(':|\.',time)
        time_converted = int(times[0])*60+int(times[1])

        if time_converted > time_last:
            quarter = quarter+1

        #print("%i, %i, %i" % ((4-quarter)*720+time_converted,time_converted,quarter))
        time_last = time_converted

        ## Find score for home and away
        score = children[3].string
        scores = score.split("-")

        writer.writerow([(4-quarter)*720+time_converted,scores[0],scores[1]])
    file.close()
