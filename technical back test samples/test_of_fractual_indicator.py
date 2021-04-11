# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 16:37:42 2021

@author: ASUS
"""

import pybbg 
import pandas as pd
from datetime import datetime, date, time
from pandas.tseries.offsets import BDay
import numpy as np
from scipy.signal import argrelextrema
from tqdm import tqdm
import matplotlib.pyplot as plt

from various_ma_package import fractal_indicator,adder

today = datetime.today()
bbg = pybbg.Pybbg()

fld_list = ['PX_OPEN', 'PX_HIGH', 'PX_LOW', 'PX_LAST']
instrument_list = ['EURUSD Curncy', 'GBPUSD Curncy', 'AUDUSD Curncy', 'USDJPY Curncy']
ticker = ['EURUSD Curncy']


trend = 10
barrier = 1
def signal(Data, what, closing, buy, sell):
    
    for i in range(len(Data)):
            
     if Data[i, what] < barrier and Data[i, closing] < Data[i - trend, closing]:
        Data[i, buy] = 1
            
     if Data[i, what] < barrier and Data[i, closing] > Data[i - trend, closing]:  
        Data[i, sell] = -1
        
    return Data
# The trend variable is the algorithm's way to see whether the market price has been trending down or up. This is to known what position to initiate (long/short) as the Fractal Indicator only shows a uniform signal which is the event of reaching 1.00
# The Data variable refers to the OHLC array
# The what variable refers to the Fractal Indicator
# The closing variable refers to the closing price
# The buy variable refers to where we should place long orders
# The sell variable refers to where we should place short orders



for ticker in instrument_list:
    #fld_list = ['open', 'close', 'high', 'low']
    #prev = (today - BDay(5))
    #my_data = bbg.bdib(ticker, fld_list, datetime(2021,1,1,1,30), datetime(2021,2,1,21,30), eventType='TRADE', interval = 30)
    my_data=bbg.bdh(ticker, fld_list, '20131201', end_date=date.today().strftime('%Y%m%d'), periodselection = 'DAILY')
    #my_data = bbg.bdib(ticker, fld_list, prev, today, eventType='TRADE', interval = 60)
    print(my_data)
    
    ewa_lookback = 20
    max_min = 14
    
    # Using the function to add 20 columns
    my_data = adder(my_data, 20)
    my_data = fractal_indicator(my_data, 1,2,ewa_lookback, max_min,4)
    my_data = signal(my_data,13,3,14,15)
    data = pd.DataFrame(my_data)
    buys = data.loc[(data[14] == 1)]
    sells = data.loc[(data[15]==-1)]
    
    plt.figure()
    plt.plot(data.index, data[3], label=ticker)
    #plt.plot(data.index, data[13], label = 'Fraction index')
    plt.plot(buys.index, data.iloc[buys.index][3], '^', markersize=10, color='g')
    plt.plot(sells.index, data.iloc[sells.index][3], 'v', markersize=10, color='r')
    plt.legend()
    plt.show()
    
    '''
    
    below is to calculate the performance
    
    
    '''
    my_data[:,16] = my_data[:,14]+my_data[:,15]
    data = pd.DataFrame(my_data)
    data[16].replace(to_replace=0, method='ffill', inplace=True)
    
    plt.figure()
    
    
    fx = data
    fx['next_day'] = data[3].shift(-1)
    fx['ret'] = fx['next_day']/data[3] - 1 
    
    (fx['ret']+1).cumprod().plot(figsize = (16,8),label='base '+ticker)
    (fx[16]*fx['ret']+1).cumprod().plot(label='strategy')
    plt.legend()
    plt.show()
