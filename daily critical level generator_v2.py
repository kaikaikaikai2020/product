# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 15:07:54 2021

@author: ASUS
"""
import pybbg 
import pandas as pd
from datetime import datetime, date, time
from pandas.tseries.offsets import BDay
import numpy as np
from scipy.signal import argrelextrema

today = datetime.today()

#prev = (today - BDay(10))
prev = (today - BDay(5))


result = pd.DataFrame(columns=['Bloomberg Code','symbol','level'])
result_low = pd.DataFrame(columns=['Bloomberg Code','symbol','level'])
result_high_v2 = pd.DataFrame(columns=['Bloomberg Code','symbol','level'])
result_low_v2 = pd.DataFrame(columns=['Bloomberg Code','symbol','level'])
from scipy.signal import find_peaks
from tqdm import tqdm
import matplotlib.pyplot as plt

def findMiddle(input_list):
    middle = float(len(input_list))/2
    if middle % 2 != 0:
        return input_list[int(middle - .5)]
    else:
        return input_list[int(middle)]

bbg = pybbg.Pybbg()

###################################################
# read TW files
###################################################
data_tw = pd.read_excel('data_tw.xlsx')
#print(data_tw)
stock_list = data_tw['Ticker']
fld_list = ['PX_OPEN', 'PX_HIGH', 'PX_LOW', 'PX_LAST']
gre='>'
les='<'
up_gap = 1
low_gap = 1
#stock_list = ['2330 TT Equity']
#data_pd = bbg.bdh(stock_list, fld_list, '20201201', end_date=date.today().strftime('%Y%m%d'), periodselection = 'DAILY')
#print(data_pd)
print('##############################################')
print('########## Start to process TW Stocks ########')
for ticker in tqdm(stock_list):
    fld_list = ['open', 'close', 'high', 'low']
    #data_pd = bbg.bdib(ticker, fld_list, datetime(2021,1,1,1,30), datetime(2021,2,1,21,30), eventType='TRADE', interval = 30)
    data_pd = bbg.bdib(ticker, fld_list, prev, today, eventType='TRADE', interval = 240)
    current = data_pd['close'].iloc[-1]*up_gap
    maxm = argrelextrema(data_pd['high'].values, np.greater)  # (array([1, 3, 6]),)
    temp_high = data_pd['high'].values[maxm]
    temp_high = temp_high[np.where(temp_high>current)]
    temp_high = np.unique(temp_high)
    temp_high = np.sort(temp_high)

    if(len(temp_high)>1):
        temp = findMiddle(list(temp_high))
        #result = result.append({'Bloomberg Code': ticker,'symbol': '>', 'level': temp_high[-1]}, ignore_index=True)
        result = result.append({'Bloomberg Code': ticker,'symbol': '>', 'level': temp}, ignore_index=True)
    
    current = data_pd['close'].iloc[-1]*low_gap
    minm = argrelextrema(data_pd['low'].values, np.less)  # (array([2, 5, 7]),)
    temp_low = data_pd['low'].values[minm]
    temp_low = temp_low[np.where(temp_low<current)]
    temp_low = np.unique(temp_low)
    temp_low = np.sort(temp_low)
    if(len(temp_low)>1):
        temp = findMiddle(list(temp_low))
        #result_low = result_low.append({'Bloomberg Code': ticker,'symbol': '<', 'level': temp_low[-1]}, ignore_index=True)    
        result_low = result_low.append({'Bloomberg Code': ticker,'symbol': '<', 'level': temp}, ignore_index=True)    

    #print('----------')
    #print('maxm',maxm)
    #print('minm',minm)
    #print(data_pd['close'].values[maxm])
    #print(data_pd['close'].values[minm])
    #print('----------')
    #peaks, _ = find_peaks(data_pd['close'].values)
    #plt.plot(x-1, y)
    #plt.plot(peaks, data_pd['close'].values[peaks], "x")
    #print(peaks)
    #print(data_pd['close'].values[peaks])
    #plt.show()
    
#print(result)

print('########## End process TW Stocks ########')


###################################################
# read KR files
###################################################
    

data_tw = pd.read_excel('data_kr.xlsx')
#print(data_tw)
stock_list = data_tw['Ticker']
fld_list = ['PX_OPEN', 'PX_HIGH', 'PX_LOW', 'PX_LAST']
#stock_list = ['2330 TT Equity']
#data_pd = bbg.bdh(stock_list, fld_list, '20201201', end_date=date.today().strftime('%Y%m%d'), periodselection = 'DAILY')
#print(data_pd)
print('##############################################')
print('########## Start to process KR Stocks ########')

for ticker in tqdm(stock_list):
    fld_list = ['open', 'close', 'high', 'low']
    #data_pd = bbg.bdib(ticker, fld_list, datetime(2021,1,1,1,30), datetime(2021,2,1,21,30), eventType='TRADE', interval = 30)
    data_pd = bbg.bdib(ticker, fld_list, prev, today, eventType='TRADE', interval = 240)
    current = data_pd['close'].iloc[-1]*up_gap
    maxm = argrelextrema(data_pd['high'].values, np.greater)  # (array([1, 3, 6]),)
    temp_high = data_pd['high'].values[maxm]
    temp_high = temp_high[np.where(temp_high>current)]
    temp_high = np.unique(temp_high)
    temp_high = np.sort(temp_high)

    if(len(temp_high)>1):
        temp = findMiddle(list(temp_high))
        #result = result.append({'Bloomberg Code': ticker,'symbol': '>', 'level': temp_high[-1]}, ignore_index=True)
        result = result.append({'Bloomberg Code': ticker,'symbol': '>', 'level': temp}, ignore_index=True)
    
    current = data_pd['close'].iloc[-1]*low_gap
    minm = argrelextrema(data_pd['low'].values, np.less)  # (array([2, 5, 7]),)
    temp_low = data_pd['low'].values[minm]
    temp_low = temp_low[np.where(temp_low<current)]
    temp_low = np.unique(temp_low)
    temp_low = np.sort(temp_low)
    if(len(temp_low)>1):
        temp = findMiddle(list(temp_low))
        #result_low = result_low.append({'Bloomberg Code': ticker,'symbol': '<', 'level': temp_low[-1]}, ignore_index=True)    
        result_low = result_low.append({'Bloomberg Code': ticker,'symbol': '<', 'level': temp}, ignore_index=True)    

    
    
    #print('----------')
    #print('maxm',maxm)
    #print('minm',minm)
    #print(data_pd['close'].values[maxm])
    #print(data_pd['close'].values[minm])
    #print('----------')
    #peaks, _ = find_peaks(data_pd['close'].values)
    #plt.plot(x-1, y)
    #plt.plot(peaks, data_pd['close'].values[peaks], "x")
    #print(peaks)
    #print(data_pd['close'].values[peaks])
    #plt.show()
print('########## End process KR Stocks ########')

###################################################
# read JP files
###################################################
    

data_tw = pd.read_excel('data_jp.xlsx')
#print(data_tw)
stock_list = data_tw['Ticker']
fld_list = ['PX_OPEN', 'PX_HIGH', 'PX_LOW', 'PX_LAST']
#stock_list = ['2330 TT Equity']
#data_pd = bbg.bdh(stock_list, fld_list, '20201201', end_date=date.today().strftime('%Y%m%d'), periodselection = 'DAILY')
#print(data_pd)
print('##############################################')
print('########## Start to process JP Stocks ########')

for ticker in tqdm(stock_list):
    fld_list = ['open', 'close', 'high', 'low']
    #data_pd = bbg.bdib(ticker, fld_list, datetime(2021,1,1,1,30), datetime(2021,2,1,21,30), eventType='TRADE', interval = 30)
    data_pd = bbg.bdib(ticker, fld_list, prev, today, eventType='TRADE', interval = 240)
    current = data_pd['close'].iloc[-1]*up_gap
    maxm = argrelextrema(data_pd['high'].values, np.greater)  # (array([1, 3, 6]),)
    temp_high = data_pd['high'].values[maxm]
    temp_high = temp_high[np.where(temp_high>current)]
    temp_high = np.unique(temp_high)
    temp_high = np.sort(temp_high)

    if(len(temp_high)>1):
        temp = findMiddle(list(temp_high))
        #result = result.append({'Bloomberg Code': ticker,'symbol': '>', 'level': temp_high[-1]}, ignore_index=True)
        result = result.append({'Bloomberg Code': ticker,'symbol': '>', 'level': temp}, ignore_index=True)
    
    current = data_pd['close'].iloc[-1]*low_gap
    minm = argrelextrema(data_pd['low'].values, np.less)  # (array([2, 5, 7]),)
    temp_low = data_pd['low'].values[minm]
    temp_low = temp_low[np.where(temp_low<current)]
    temp_low = np.unique(temp_low)
    temp_low = np.sort(temp_low)
    if(len(temp_low)>1):
        temp = findMiddle(list(temp_low))
        #result_low = result_low.append({'Bloomberg Code': ticker,'symbol': '<', 'level': temp_low[-1]}, ignore_index=True)    
        result_low = result_low.append({'Bloomberg Code': ticker,'symbol': '<', 'level': temp}, ignore_index=True)    

    
    
    #print('----------')
    #print('maxm',maxm)
    #print('minm',minm)
    #print(data_pd['close'].values[maxm])
    #print(data_pd['close'].values[minm])
    #print('----------')
    #peaks, _ = find_peaks(data_pd['close'].values)
    #plt.plot(x-1, y)
    #plt.plot(peaks, data_pd['close'].values[peaks], "x")
    #print(peaks)
    #print(data_pd['close'].values[peaks])
    #plt.show()

print('########## End process JP Stocks ########')
###################################################
# read HK files
###################################################
    

data_tw = pd.read_excel('data_hk.xlsx')
#print(data_tw)
stock_list = data_tw['Ticker']
fld_list = ['PX_OPEN', 'PX_HIGH', 'PX_LOW', 'PX_LAST']
#stock_list = ['2330 TT Equity']
#data_pd = bbg.bdh(stock_list, fld_list, '20201201', end_date=date.today().strftime('%Y%m%d'), periodselection = 'DAILY')
#print(data_pd)
print('##############################################')
print('########## Start to process HK Stocks ########')

for ticker in tqdm(stock_list):
    fld_list = ['open', 'close', 'high', 'low']
    #data_pd = bbg.bdib(ticker, fld_list, datetime(2021,1,1,1,30), datetime(2021,2,1,21,30), eventType='TRADE', interval = 30)
    data_pd = bbg.bdib(ticker, fld_list, prev, today, eventType='TRADE', interval = 240)
    current = data_pd['close'].iloc[-1]*up_gap
    maxm = argrelextrema(data_pd['high'].values, np.greater)  # (array([1, 3, 6]),)
    temp_high = data_pd['high'].values[maxm]
    temp_high = temp_high[np.where(temp_high>current)]
    temp_high = np.unique(temp_high)
    temp_high = np.sort(temp_high)

    if(len(temp_high)>1):
        temp = findMiddle(list(temp_high))
        #result = result.append({'Bloomberg Code': ticker,'symbol': '>', 'level': temp_high[-1]}, ignore_index=True)
        result = result.append({'Bloomberg Code': ticker,'symbol': '>', 'level': temp}, ignore_index=True)
    
    current = data_pd['close'].iloc[-1]*low_gap
    minm = argrelextrema(data_pd['low'].values, np.less)  # (array([2, 5, 7]),)
    temp_low = data_pd['low'].values[minm]
    temp_low = temp_low[np.where(temp_low<current)]
    temp_low = np.unique(temp_low)
    temp_low = np.sort(temp_low)
    if(len(temp_low)>1):
        temp = findMiddle(list(temp_low))
        #result_low = result_low.append({'Bloomberg Code': ticker,'symbol': '<', 'level': temp_low[-1]}, ignore_index=True)    
        result_low = result_low.append({'Bloomberg Code': ticker,'symbol': '<', 'level': temp}, ignore_index=True)    

    
    
    #print('----------')
    #print('maxm',maxm)
    #print('minm',minm)
    #print(data_pd['close'].values[maxm])
    #print(data_pd['close'].values[minm])
    #print('----------')
    #peaks, _ = find_peaks(data_pd['close'].values)
    #plt.plot(x-1, y)
    #plt.plot(peaks, data_pd['close'].values[peaks], "x")
    #print(peaks)
    #print(data_pd['close'].values[peaks])
    #plt.show()

print('########## End process HK Stocks ########')

###################################################
# read SZ files
###################################################
    

data_tw = pd.read_excel('data_sz.xlsx')
#print(data_tw)
stock_list = data_tw['Ticker']
fld_list = ['PX_OPEN', 'PX_HIGH', 'PX_LOW', 'PX_LAST']
#stock_list = ['2330 TT Equity']
#data_pd = bbg.bdh(stock_list, fld_list, '20201201', end_date=date.today().strftime('%Y%m%d'), periodselection = 'DAILY')
#print(data_pd)
print('##############################################')
print('########## Start to process SZ Stocks ########')

for ticker in tqdm(stock_list):
    fld_list = ['open', 'close', 'high', 'low']
    #data_pd = bbg.bdib(ticker, fld_list, datetime(2021,1,1,1,30), datetime(2021,2,1,21,30), eventType='TRADE', interval = 30)
    data_pd = bbg.bdib(ticker, fld_list, prev, today, eventType='TRADE', interval = 240)
    current = data_pd['close'].iloc[-1]*up_gap
    maxm = argrelextrema(data_pd['high'].values, np.greater)  # (array([1, 3, 6]),)
    temp_high = data_pd['high'].values[maxm]
    temp_high = temp_high[np.where(temp_high>current)]
    temp_high = np.unique(temp_high)
    temp_high = np.sort(temp_high)

    if(len(temp_high)>1):
        temp = findMiddle(list(temp_high))
        #result = result.append({'Bloomberg Code': ticker,'symbol': '>', 'level': temp_high[-1]}, ignore_index=True)
        result = result.append({'Bloomberg Code': ticker,'symbol': '>', 'level': temp}, ignore_index=True)
    
    current = data_pd['close'].iloc[-1]*low_gap
    minm = argrelextrema(data_pd['low'].values, np.less)  # (array([2, 5, 7]),)
    temp_low = data_pd['low'].values[minm]
    temp_low = temp_low[np.where(temp_low<current)]
    temp_low = np.unique(temp_low)
    temp_low = np.sort(temp_low)
    if(len(temp_low)>1):
        temp = findMiddle(list(temp_low))
        #result_low = result_low.append({'Bloomberg Code': ticker,'symbol': '<', 'level': temp_low[-1]}, ignore_index=True)    
        result_low = result_low.append({'Bloomberg Code': ticker,'symbol': '<', 'level': temp}, ignore_index=True)    
    
    #print('----------')
    #print('maxm',maxm)
    #print('minm',minm)
    #print(data_pd['close'].values[maxm])
    #print(data_pd['close'].values[minm])
    #print('----------')
    #peaks, _ = find_peaks(data_pd['close'].values)
    #plt.plot(x-1, y)
    #plt.plot(peaks, data_pd['close'].values[peaks], "x")
    #print(peaks)
    #print(data_pd['close'].values[peaks])
    #plt.show()
print('########## End process SZ Stocks ########')
###################################################
# read SH files
###################################################
    

data_tw = pd.read_excel('data_sh.xlsx')
#print(data_tw)
stock_list = data_tw['Ticker']
fld_list = ['PX_OPEN', 'PX_HIGH', 'PX_LOW', 'PX_LAST']
#stock_list = ['2330 TT Equity']
#data_pd = bbg.bdh(stock_list, fld_list, '20201201', end_date=date.today().strftime('%Y%m%d'), periodselection = 'DAILY')
#print(data_pd)
print('##############################################')
print('########## Start to process SH Stocks ########')
for ticker in tqdm(stock_list):
    fld_list = ['open', 'close', 'high', 'low']
    #data_pd = bbg.bdib(ticker, fld_list, datetime(2021,1,1,1,30), datetime(2021,2,1,21,30), eventType='TRADE', interval = 30)
    data_pd = bbg.bdib(ticker, fld_list, prev, today, eventType='TRADE', interval = 240)
    current = data_pd['close'].iloc[-1]*up_gap
    maxm = argrelextrema(data_pd['high'].values, np.greater)  # (array([1, 3, 6]),)
    temp_high = data_pd['high'].values[maxm]
    temp_high = temp_high[np.where(temp_high>current)]
    temp_high = np.unique(temp_high)
    temp_high = np.sort(temp_high)

    if(len(temp_high)>1):
        temp = findMiddle(list(temp_high))
        #result = result.append({'Bloomberg Code': ticker,'symbol': '>', 'level': temp_high[-1]}, ignore_index=True)
        result = result.append({'Bloomberg Code': ticker,'symbol': '>', 'level': temp}, ignore_index=True)
    
    current = data_pd['close'].iloc[-1]*low_gap
    minm = argrelextrema(data_pd['low'].values, np.less)  # (array([2, 5, 7]),)
    temp_low = data_pd['low'].values[minm]
    temp_low = temp_low[np.where(temp_low<current)]
    temp_low = np.unique(temp_low)
    temp_low = np.sort(temp_low)
    if(len(temp_low)>1):
        temp = findMiddle(list(temp_low))
        #result_low = result_low.append({'Bloomberg Code': ticker,'symbol': '<', 'level': temp_low[-1]}, ignore_index=True)    
        result_low = result_low.append({'Bloomberg Code': ticker,'symbol': '<', 'level': temp}, ignore_index=True)    

    
    #print('----------')
    #print('maxm',maxm)
    #print('minm',minm)
    #print(data_pd['close'].values[maxm])
    #print(data_pd['close'].values[minm])
    #print('----------')
    #peaks, _ = find_peaks(data_pd['close'].values)
    #plt.plot(x-1, y)
    #plt.plot(peaks, data_pd['close'].values[peaks], "x")
    #print(peaks)
    #print(data_pd['close'].values[peaks])
    #plt.show()
print('########## End process SH Stocks ########')  


data_tw = pd.read_excel('data_us.xlsx')
#print(data_tw)
stock_list = data_tw['Ticker']
fld_list = ['PX_OPEN', 'PX_HIGH', 'PX_LOW', 'PX_LAST']


print('##############################################')
print('########## Start to process US Stocks ########')
for ticker in tqdm(stock_list):
    fld_list = ['open', 'close', 'high', 'low']
    #data_pd = bbg.bdib(ticker, fld_list, datetime(2021,1,1,1,30), datetime(2021,2,1,21,30), eventType='TRADE', interval = 30)
    data_pd = bbg.bdib(ticker, fld_list, prev, today, eventType='TRADE', interval = 240)
    current = data_pd['close'].iloc[-1]*up_gap
    maxm = argrelextrema(data_pd['high'].values, np.greater)  # (array([1, 3, 6]),)
    temp_high = data_pd['high'].values[maxm]
    temp_high = temp_high[np.where(temp_high>current)]
    temp_high = np.unique(temp_high)
    temp_high = np.sort(temp_high)

    if(len(temp_high)>1):
        temp = findMiddle(list(temp_high))
        #result = result.append({'Bloomberg Code': ticker,'symbol': '>', 'level': temp_high[-1]}, ignore_index=True)
        result = result.append({'Bloomberg Code': ticker,'symbol': '>', 'level': temp}, ignore_index=True)
    
    current = data_pd['close'].iloc[-1]*low_gap
    minm = argrelextrema(data_pd['low'].values, np.less)  # (array([2, 5, 7]),)
    temp_low = data_pd['low'].values[minm]
    temp_low = temp_low[np.where(temp_low<current)]
    temp_low = np.unique(temp_low)
    temp_low = np.sort(temp_low)
    if(len(temp_low)>1):
        temp = findMiddle(list(temp_low))
        #result_low = result_low.append({'Bloomberg Code': ticker,'symbol': '<', 'level': temp_low[-1]}, ignore_index=True)    
        result_low = result_low.append({'Bloomberg Code': ticker,'symbol': '<', 'level': temp}, ignore_index=True)    

    
    #print('----------')
    #print('maxm',maxm)
    #print('minm',minm)
    #print(data_pd['close'].values[maxm])
    #print(data_pd['close'].values[minm])
    #print('----------')
    #peaks, _ = find_peaks(data_pd['close'].values)
    #plt.plot(x-1, y)
    #plt.plot(peaks, data_pd['close'].values[peaks], "x")
    #print(peaks)
    #print(data_pd['close'].values[peaks])
    #plt.show()
print('########## End process SH Stocks ########')  



today = today.strftime('%Y-%m-%d')
filename = 'all_stock {0}_high.csv'.format(today)  
result.to_csv(filename)
filename = 'all_stock {0}_low.csv'.format(today)  
result_low.to_csv(filename)

print(result)

    
    
