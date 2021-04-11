# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 09:24:36 2021

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

today = datetime.today()
bbg = pybbg.Pybbg()

fld_list = ['PX_OPEN', 'PX_HIGH', 'PX_LOW', 'PX_LAST']
instrument_list = ['EURUSD Curncy', 'GBPUSD Curcny', 'AUDUSD Curncy', 'USDJPY Curncy']
ticker = ['EURUSD Curncy']

#fld_list = ['open', 'close', 'high', 'low']
#prev = (today - BDay(5))
#my_data = bbg.bdib(ticker, fld_list, datetime(2021,1,1,1,30), datetime(2021,2,1,21,30), eventType='TRADE', interval = 30)
my_data=bbg.bdh(ticker, fld_list, '20131201', end_date=date.today().strftime('%Y%m%d'), periodselection = 'DAILY')
#my_data = bbg.bdib(ticker, fld_list, prev, today, eventType='TRADE', interval = 60)
print(my_data)

def donchian(Data, low, high, where_up, where_down, median, period):
    
    for i in range(len(Data)):
        try:
            Data[i, where_up] = max(Data[i - period:i + 1, 1])
        except ValueError:
            pass
    for i in range(len(Data)):
        try:
            Data[i, where_down] = min(Data[i - period:i + 1, 2]) 
        except ValueError:
            pass        
    for i in range(len(Data)): 
        try:
            Data[i, median] = (Data[i, where_up] + Data[i, where_down]) / 2 
        except ValueError:
            pass
    return Data


def adder(Data, times):
    for i in range(1, times + 1):
        z = np.zeros((len(Data), 1), dtype = float)
        Data = np.append(Data, z, axis = 1)
    return Data
# Using the function to add 20 columns
my_data = adder(my_data, 20)
my_data = donchian(my_data, 2, 1, 4, 5, 6, 20)

'''
plt.plot(my_data[-500:, 3], color = 'black')
plt.plot(my_data[-500:, 4])
plt.plot(my_data[-500:, 5])
plt.plot(my_data[-500:, 6])
plt.grid()
'''
where_close = 3

def donchian_signals(Data, onwhat, donch_up, donch_down):
    
    for i in range(len(Data)):
        
          if Data[i, where_close] > Data[i - 1, donch_up]:
              Data[i, onwhat] = 1
         
          elif Data[i, where_close] < Data[i - 1, donch_down]:
              Data[i, onwhat+1] = -1 
              
    return Data
# Note: The onwhat variable is the closing price
my_data = donchian_signals(my_data,7,4,5)


data = pd.DataFrame(my_data)
buys = data.loc[(data[7] == 1)]
sells = data.loc[(data[8]==-1)]

plt.plot(data.index, data[3], label='Euro')
plt.plot(buys.index, data.iloc[buys.index][3], '^', markersize=10, color='g')
plt.plot(sells.index, data.iloc[sells.index][3], 'v', markersize=10, color='r')




