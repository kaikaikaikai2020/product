# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 23:07:19 2021

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

from various_ma_package import volatility,adder

today = datetime.today()
bbg = pybbg.Pybbg()

fld_list = ['PX_OPEN', 'PX_HIGH', 'PX_LOW', 'PX_LAST']
instrument_list = ['EURUSD Curncy', 'GBPUSD Curncy', 'AUDUSD Curncy', 'USDJPY Curncy']
ticker = ['EURUSD Curncy']



def vri(Data, lookback, what, high, low, where):
    # Volatility Calculation
    Data = volatility(Data, lookback, what, where)
    # Differencing as a Measure of Momentum
    for i in range(len(Data)):    
        Data[i, where + 1] = Data[i, what] - Data[i - lookback, 0]
    # Extreme Range
    for i in range(len(Data)):
        try:
            Data[i, where + 2] = max(Data[i - lookback:i + 1, high])
        except ValueError:
            pass
    for i in range(len(Data)):
        try:
            Data[i, where + 3] = min(Data[i - lookback:i + 1, low])
        except ValueError:
            pass
        
    Data[:, where + 4] = Data[:, where + 1] / (Data[:, where + 2] - Data[:, where + 3])    
        
    Data[:, where + 4] = Data[:, where + 4] * Data[:, where] * 1000
       
    return Data

def signal(Data, what, barrier, buy, sell):
    
    for i in range(len(Data)):
            
     if Data[i, what] > -barrier and Data[i-1, what] < -barrier and Data[i-2, what]<-barrier:
        Data[i, buy] = 1
            
     if Data[i, what] < barrier and Data[i-1, what] > barrier and Data[i-2, what] > barrier:  
        Data[i, sell] = -1
        
    return Data

for ticker in instrument_list:
    my_data=bbg.bdh(ticker, fld_list, '20131201', end_date=date.today().strftime('%Y%m%d'), periodselection = 'DAILY')
    print(my_data)
    
    
    # Using the function to add 20 columns
    my_data = adder(my_data, 20)
    my_data = vri(my_data, 3,3,1,2,4)
    my_data = signal(my_data,8,3,9,10)
    data = pd.DataFrame(my_data)
    #data[8].plot()
    buys = data.loc[(data[9] == 1)]
    sells = data.loc[(data[10]==-1)]
    
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
    my_data[:,11] = my_data[:,9]+my_data[:,10]
    data = pd.DataFrame(my_data)
    data[11].replace(to_replace=0, method='ffill', inplace=True)
    
    plt.figure()
    
    
    fx = data
    fx['next_day'] = data[3].shift(-1)
    fx['ret'] = fx['next_day']/data[3] - 1 
    
    (fx['ret']+1).cumprod().plot(figsize = (16,8),label='base '+ticker)
    (fx[11]*fx['ret']+1).cumprod().plot(label='strategy')
    plt.legend()
    plt.show()
