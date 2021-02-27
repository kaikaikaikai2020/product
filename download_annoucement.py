# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 19:48:59 2020

@author: ASUS
"""

import uqer
from uqer import DataAPI
import pandas as pd
client = uqer.Client(token='473c6c98362ce267eafb83a9bf37bcd64eb25d739829b1356c426bc1d3f35ac7')

stock_code = DataAPI.IdxConsGet(secID=u"",ticker=u"399106",isNew=u"",intoDate=u"20201218",field=u"",pandas="1")
#print(stock_code)

df = pd.DataFrame()
for code in stock_code['consTickerSymbol']:
    print(code)
    df = df.append(DataAPI.AnnoProfitabilityGet(ticker=code,reportID=u"",beginDate=u"20100101",endDate=u"20201218",field=u"",pandas="1"))
    #print(df)
    
df.to_csv('shenzhen announcement.csv')