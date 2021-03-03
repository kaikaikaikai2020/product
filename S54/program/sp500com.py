# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 22:36:07 2020

@author: adair2019
"""

import bs4 as bs
import pickle
import requests
def save_sp500_tickers():
   resp = requests.get('http://en.wikipedia.org/wiki/List_of_S%26P_500_companies')
   soup = bs.BeautifulSoup(resp.text, 'lxml')
   table = soup.find('table', {'class': 'wikitable sortable'})
   tickers = []
   for row in table.findAll('tr')[1:]:
       ticker = row.findAll('td')[0].text
       tickers.append(ticker)
       
   with open("sp500tickers.pickle","wb") as f:
       pickle.dump(tickers,f)
       
   return tickers

save_sp500_tickers()# 获取S&P500 里面的内容 zhujin 注
with open('sp500tickers.pickle', 'rb') as f:
   # The protocol version used is detected automatically, so we do not
   # have to specify it.
   data = pickle.load(f)
   
print('Top 100 list of S&P500 %s' %data[:100])