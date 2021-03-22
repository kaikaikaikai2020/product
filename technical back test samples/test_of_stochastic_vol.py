# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 14:29:22 2021

@author: ASUS
"""
from various_ma_package import volatility,normalizer,jump,adder

import pybbg 
import pandas as pd
from datetime import datetime, date, time
from pandas.tseries.offsets import BDay
import numpy as np
from scipy.signal import argrelextrema
from tqdm import tqdm
import matplotlib.pyplot as plt



def fib_stoch(Data, volatility_lookback, what, where):
        
    Data = volatility(Data, volatility_lookback, what, where)
    Data = normalizer(Data, volatility_lookback, where, where + 1)
    for i in range(len(Data)):
        Data[i, where + 1] = round(Data[i, where + 1], 0) 
    
    for i in range(len(Data)):
        
        if Data[i, where + 1] >= 0 and Data[i, where + 1] <= 10 :
            Data[i, where + 1] = 144
        if Data[i, where + 1] > 10 and Data[i, where + 1] <= 20 :
            Data[i, where + 1] = 89            
        if Data[i, where + 1] > 20 and Data[i, where + 1] <= 30 :
            Data[i, where + 1] = 55
        if Data[i, where + 1] > 30 and Data[i, where + 1] <= 40 :
            Data[i, where + 1] = 34
        if Data[i, where + 1] > 40 and Data[i, where + 1] <= 50 :
            Data[i, where + 1] = 21
        if Data[i, where + 1] > 50 and Data[i, where + 1] <= 60 :
            Data[i, where + 1] = 13
        if Data[i, where + 1] > 60 and Data[i, where + 1] <= 70 :
            Data[i, where + 1] = 8
        if Data[i, where + 1] > 70 and Data[i, where + 1] <= 80 :
            Data[i, where + 1] = 5
        if Data[i, where + 1] > 80 and Data[i, where + 1] <= 90 :
            Data[i, where + 1] = 3
        if Data[i, where + 1] > 90 and Data[i, where + 1] <= 100 :
            Data[i, where + 1] = 2            
    
    Data = jump(Data, volatility_lookback)
        
    for i in range(len(Data)):
        
        try:
            lookback = int(Data[i, where + 1])
            Data[i, where + 2] = (Data[i, what] - min(Data[i - lookback + 1:i + 1, 2])) / (max(Data[i - lookback + 1:i + 1, 1]) -    min(Data[i - lookback + 1:i + 1, 2]))
        
        except ValueError:
            pass
    return Data


'''

Go long (Buy) whenever the oscillator reaches 10 with the previous value being above 10. Hold this position until getting another signal or getting closed out due to risk management mentioned above.
Go short (Sell) whenever the oscillator reaches 90 with the previous value being less than 90. Hold this position until getting another signal or getting closed out due to risk management mentioned above.


'''
today = datetime.today()
bbg = pybbg.Pybbg()

fld_list = ['PX_OPEN', 'PX_HIGH', 'PX_LOW', 'PX_LAST']
instrument_list = ['EURUSD Curncy', 'GBPUSD Curncy', 'AUDUSD Curncy', 'USDJPY Curncy']
ticker = ['EURUSD Curncy']


for ticker in instrument_list:
    #fld_list = ['open', 'close', 'high', 'low']
    #prev = (today - BDay(5))
    #my_data = bbg.bdib(ticker, fld_list, datetime(2021,1,1,1,30), datetime(2021,2,1,21,30), eventType='TRADE', interval = 30)
    my_data=bbg.bdh(ticker, fld_list, '20131201', end_date=date.today().strftime('%Y%m%d'), periodselection = 'DAILY')
    #my_data = bbg.bdib(ticker, fld_list, prev, today, eventType='TRADE', interval = 60)
    print(my_data)
    
    lookback = 5
    up_barrier = 90
    low_barrier = 10
    width = 20
    
    # Using the function to add 20 columns
    my_data = adder(my_data, 20)
    my_data = fib_stoch(my_data, lookback,3,4)
    data = pd.DataFrame(my_data)
    buys = data.loc[(data[6] == 1)]
    sells = data.loc[(data[7]==-1)]
    
    plt.figure()
    plt.plot(data.index, data[3], label='Euro')
    plt.plot(buys.index, data.iloc[buys.index][3], '^', markersize=10, color='g')
    plt.plot(sells.index, data.iloc[sells.index][3], 'v', markersize=10, color='r')
    
    '''
    
    below is to calculate the performance
    
    
    '''
    my_data[:,8] = my_data[:,6]+my_data[:,7]
    data = pd.DataFrame(my_data)
    data[8].replace(to_replace=0, method='ffill', inplace=True)
    
    plt.figure()
    
    
    fx = data
    fx['next_day'] = data[3].shift(-1)
    fx['ret'] = fx['next_day']/data[3] - 1 
    
    (fx['ret']+1).cumprod().plot(figsize = (16,8),label='base '+ticker)
    (fx[8]*fx['ret']+1).cumprod().plot(label='strategy')
    plt.legend()
    plt.show()
