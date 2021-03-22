# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 13:55:22 2021

@author: ASUS
"""
from various_ma_package import rolling_correlation,rsi, adder

import pybbg 
import pandas as pd
from datetime import datetime, date, time
from pandas.tseries.offsets import BDay
import numpy as np
from scipy.signal import argrelextrema
from tqdm import tqdm
import matplotlib.pyplot as plt



today = datetime.today()
bbg = pybbg.Pybbg()

fld_list = ['PX_OPEN', 'PX_HIGH', 'PX_LOW', 'PX_LAST']
instrument_list = ['EURUSD Curncy', 'GBPUSD Curncy', 'AUDUSD Curncy', 'USDJPY Curncy']
ticker = ['EURUSD Curncy']


momentum_lookback = 14
corr_lookback =10
up_barrier = 70
low_barrier = 30


def dynamic_rsi(Data, momentum_lookback, corr_lookback, what, where):
        
    for i in range(len(Data)):
        Data[i, where] = Data[i, what] / Data[i - momentum_lookback, what] * 100
    Data = rolling_correlation(Data, what, where, corr_lookback, where + 1)
    
    for i in range(len(Data)):
        
        if Data[i, where + 1]>= -1.00 and Data[i, where + 1]<= 0.10:
            Data[i, where + 1] = 14
        if Data[i, where + 1] > 0.10 and Data[i, where + 1]<= 0.20:
            Data[i, where + 1] = 10            
        if Data[i, where + 1] > 0.20 and Data[i, where + 1]<= 0.30:
            Data[i, where + 1] = 9
        if Data[i, where + 1] > 0.30 and Data[i, where + 1]<= 0.40:
            Data[i, where + 1] = 8
        if Data[i, where + 1] > 0.40 and Data[i, where + 1]<= 0.50:
            Data[i, where + 1] = 7
        if Data[i, where + 1] > 0.50 and Data[i, where + 1]<= 0.60:
            Data[i, where + 1] = 6
        if Data[i, where + 1] > 0.60 and Data[i, where + 1]<= 0.70:
            Data[i, where + 1] = 5
        if Data[i, where + 1] > 0.70 and Data[i, where + 1]<= 0.80:
            Data[i, where + 1] = 4
        if Data[i, where + 1] > 0.80 and Data[i, where + 1]<= 0.90:
            Data[i, where + 1] = 3
        if Data[i, where + 1] > 0.90 and Data[i, where + 1]<= 1.00:
            Data[i, where + 1] = 2            
    
    Data = rsi(Data, 14, 3, 0)
    Data = rsi(Data, 10, 3, 0)
    Data = rsi(Data, 9, 3, 0)
    Data = rsi(Data, 8, 3, 0)
    Data = rsi(Data, 7, 3, 0)
    Data = rsi(Data, 6, 3, 0)
    Data = rsi(Data, 5, 3, 0)
    Data = rsi(Data, 4, 3, 0)
    Data = rsi(Data, 3, 3, 0)
    Data = rsi(Data, 2, 3, 0)
    for i in range(len(Data)):
        
        if Data[i, where + 1] == 14:
            Data[i, where + 12] = Data[i, where + 2]
        if Data[i, where + 1] == 10:
            Data[i, where + 12] = Data[i, where + 3]
        if Data[i, where + 1] == 9:
            Data[i, where + 12] = Data[i, where + 4]
        if Data[i, where + 1] == 8:
            Data[i, where + 12] = Data[i, where + 5]
        if Data[i, where + 1] == 7:
            Data[i, where + 12] = Data[i, where + 6]
        if Data[i, where + 1] == 6:
            Data[i, where + 12] = Data[i, where + 7]
        if Data[i, where + 1] == 5:
            Data[i, where + 12] = Data[i, where + 8]
        if Data[i, where + 1] == 4:
            Data[i, where + 12] = Data[i, where + 9]
        if Data[i, where + 1] == 3:
            Data[i, where + 12] = Data[i, where + 10]
        if Data[i, where + 1] == 2:            
            Data[i, where + 12] = Data[i, where + 11]
    return Data

'''

The below graphs show the performance tables and the equity curves of a strategy that:
Goes long (Buy) whenever the indicator reaches 30 and hold the position until a new signal is generated or until the position is stopped-out either by a stop-loss or a profit order.
Goes short (Sell) whenever the indicator reaches 70 and hold the position until a new signal is generated or until the position is stopped-out either by a stop-loss or a profit order.


'''
def signal(Data, what, up_barrier, down_barrier, buy, sell):
    
    for i in range(len(Data)):
            
     if Data[i, what] > up_barrier:
        Data[i, buy] = -1
            
     if Data[i, what] < down_barrier:  
        Data[i, sell] = 1
        
    return Data

for ticker in instrument_list:
    #fld_list = ['open', 'close', 'high', 'low']
    #prev = (today - BDay(5))
    #my_data = bbg.bdib(ticker, fld_list, datetime(2021,1,1,1,30), datetime(2021,2,1,21,30), eventType='TRADE', interval = 30)
    my_data=bbg.bdh(ticker, fld_list, '20131201', end_date=date.today().strftime('%Y%m%d'), periodselection = 'DAILY')
    #my_data = bbg.bdib(ticker, fld_list, prev, today, eventType='TRADE', interval = 60)
    print(my_data)
    
    
    # Using the function to add 20 columns
    my_data = adder(my_data, 20)
    my_data = dynamic_rsi(my_data, momentum_lookback, corr_lookback, 3, 4)
    my_data = signal(my_data,4, up_barrier, low_barrier,5,6)
    data = pd.DataFrame(my_data)
    buys = data.loc[(data[5] == 1)]
    sells = data.loc[(data[6]==-1)]
    
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

