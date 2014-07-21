from bs4 import BeautifulSoup
from selenium import webdriver
import requests
import csv
import re
import time

driver = webdriver.Firefox()
for league in ["uefaeuropaleague","uefachampionsleague"]:
    for year in range(1994,2014):
        print(year)
        driver.get('http://www.uefa.com/%s/season=%s/matches/all/index.html' % (league,year))
        time.sleep(10)
        html = driver.page_source
        #soup = BeautifulSoup(open("/Users/afranks/Dropbox/tvranking/Data/UEFA/%s.html" % (year),encoding="utf8"))
        soup = BeautifulSoup(html)
        dates = soup.findAll('table',{'class':'tdate'})
        file = open("/Users/afranks/Dropbox/tvranking/Data/UEFA/%s-%s.csv" % (league,year),'w',encoding="utf8")
        writer = csv.writer(file)
        writer.writerow(["Year","Month","Day","Home","Away","Home Goals","Away Goals"])
        for dateTag in dates:
            date = dateTag['class'][1]
            month = date[8:10]
            day = date[10:12]
            games = dateTag.findAll('tbody')
            for game in games:
                game.find('tr',{'class':'match_res'}).find('td',{'class':'home'})
                home = game.find('tr',{'class':'match_res'}).find('td',{'class':'home'}).find('a').string
                away = game.find('tr',{'class':'match_res'}).find('td',{'class':'away'}).find('a').string
                score = game.find('tr',{'class':'match_res'}).find('td',{'class':'score'}).find('a')
                if score==None:
                    score = game.find('tr',{'class':'match_res'}).find('td',{'class':'score'})
                score = score.text
                scores = re.split('-| ',score)
                writer.writerow([year,month,day,home,away,scores[0],scores[1]])
    file.close()
driver.close()
