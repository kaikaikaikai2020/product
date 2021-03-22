# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 09:01:48 2021

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

from various_ma_package import adx,adder


today = datetime.today()
bbg = pybbg.Pybbg()

fld_list = ['PX_OPEN', 'PX_HIGH', 'PX_LOW', 'PX_LAST']
instrument_list = ['EURUSD Curncy', 'GBPUSD Curncy', 'AUDUSD Curncy', 'USDJPY Curncy']
ticker = ['EURUSD Curncy']



def signal(Data, di_plus, di_minus, buy, sell):
    
    for i in range(len(Data)):
        if Data[i, di_plus] > Data[i, di_minus] and Data[i - 1, di_plus] < Data[i - 1, di_minus]:
            Data[i, buy] = 1
            
        if Data[i, di_minus] > Data[i, di_plus] and Data[i - 1, di_minus] < Data[i - 1, di_plus]:
            Data[i, sell] = -1
            
    return Data
# The Data variable is the OHLC array
# The di_plus is the column of the +DI
# The di_minus is the column of the -DI
# The buy variable is the column where buy orders are put
# The sell variable is the column where sell orders are put


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
    my_data = adx(my_data, 1,2,3, ewa_lookback, 4)
    my_data = signal(my_data,9,10,15,16)
    data = pd.DataFrame(my_data)
    buys = data.loc[(data[15] == 1)]
    sells = data.loc[(data[16]==-1)]
    
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
    my_data[:,17] = my_data[:,15]+my_data[:,16]
    data = pd.DataFrame(my_data)
    data[17].replace(to_replace=0, method='ffill', inplace=True)
    
    plt.figure()
    
    
    fx = data
    fx['next_day'] = data[3].shift(-1)
    fx['ret'] = fx['next_day']/data[3] - 1 
    
    (fx['ret']+1).cumprod().plot(figsize = (16,8),label='base '+ticker)
    (fx[17]*fx['ret']+1).cumprod().plot(label='strategy')
    plt.legend()
    plt.show()

lower_barrier = 25
upper_barrier = 75
adx_limit = 25


def signal2(Data, rsi_column, adx_column, buy, sell):
    
    for i in range(len(Data)):
            
        if Data[i, rsi_column] < lower_barrier and Data[i - 1, rsi_column] > lower_barrier and Data[i, adx_column] < adx_limit:
            Data[i, buy] = 1
            
        if Data[i, rsi_column] > upper_barrier and Data[i - 1, rsi_column] < upper_barrier and Data[i, adx_column] < adx_limit:
            Data[i, sell] = -1
# adx_limit is a variable that you can set to represent the ADX low value


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
    my_data = adx(my_data, 1,2,3, ewa_lookback, 4)
    my_data = rsi(my_data)
    my_data = signal(my_data,9,10,15,16)
    data = pd.DataFrame(my_data)
    buys = data.loc[(data[15] == 1)]
    sells = data.loc[(data[16]==-1)]
    
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
    my_data[:,17] = my_data[:,15]+my_data[:,16]
    data = pd.DataFrame(my_data)
    data[17].replace(to_replace=0, method='ffill', inplace=True)
    
    plt.figure()
    
    
    fx = data
    fx['next_day'] = data[3].shift(-1)
    fx['ret'] = fx['next_day']/data[3] - 1 
    
    (fx['ret']+1).cumprod().plot(figsize = (16,8),label='base '+ticker)
    (fx[17]*fx['ret']+1).cumprod().plot(label='strategy')
    plt.legend()
    plt.show()
