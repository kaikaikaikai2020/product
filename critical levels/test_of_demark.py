# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 11:19:07 2021

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



def TD_differential(Data, true_low, true_high, buy, sell):
    for i in range(len(Data)):
        
        # True low
        #Data[i, true_low] = min(Data[i, 2], Data[i - 1, 3])
        Data[i, true_low] = min(Data[i, 2], Data[i - 1, 2])
        Data[i, true_low] = Data[i, 3] - Data[i, true_low]
            
        # True high  
        Data[i, true_high] = max(Data[i, 1], Data[i - 1, 1])
        Data[i, true_high] = Data[i, 3] - Data[i, true_high]
        
        # TD Differential
        if Data[i, 3] < Data[i - 1, 3] and Data[i - 1, 3] < Data[i - 2, 3] and Data[i, true_low] > Data[i - 1, true_low] and Data[i, true_high] < Data[i - 1, true_high]: 
               Data[i, buy] = 1
        if Data[i, 3] > Data[i - 1, 3] and Data[i - 1, 3] > Data[i - 2, 3] and Data[i, true_low] < Data[i - 1, true_low] and Data[i, true_high] > Data[i - 1, true_high]: 
               Data[i, sell] = -1
    return Data


               
def TD_reverse_differential(Data, true_low, true_high, buy, sell):
    for i in range(len(Data)):
        
        # True low
        Data[i, true_low] = min(Data[i, 2], Data[i - 1, 2])
        Data[i, true_low] = Data[i, 3] - Data[i, true_low]
            
        # True high  
        Data[i, true_high] = max(Data[i, 1], Data[i - 1, 1])
        Data[i, true_high] = Data[i, 3] - Data[i, true_high]
        
        # TD Differential
        if Data[i, 3] < Data[i - 1, 3] and Data[i - 1, 3] < Data[i - 2, 3] and Data[i, true_low] < Data[i - 1, true_low] and Data[i, true_high] > Data[i - 1, true_high]: 
               Data[i, buy] = 1
        if Data[i, 3] > Data[i - 1, 3] and Data[i - 1, 3] > Data[i - 2, 3] and Data[i, true_low] > Data[i - 1, true_low] and Data[i, true_high] < Data[i - 1, true_high]: 
               Data[i, sell] = -1
    return Data
           

def TD_anti_differential(Data, true_low, true_high, buy, sell):
    for i in range(len(Data)):
        
        if Data[i, 3] < Data[i - 1, 3] and Data[i - 1, 3] > Data[i - 2, 3] and Data[i - 2, 3] < Data[i - 3, 3] and Data[i - 3, 3] < Data[i - 4, 3]: 
               Data[i, buy] = 1
        if Data[i, 3] > Data[i - 1, 3] and Data[i - 1, 3] < Data[i - 2, 3] and Data[i - 2, 3] > Data[i - 3, 3] and Data[i - 3, 3] > Data[i - 4, 3]: 
               Data[i, sell] = -1               
               
    return Data

def adder(Data, times):
    for i in range(1, times + 1):
        z = np.zeros((len(Data), 1), dtype = float)
        Data = np.append(Data, z, axis = 1)
    return Data
# Using the function to add 20 columns
my_data = adder(my_data, 20)
my_data = TD_differential(my_data,4,5,6,7)


data = pd.DataFrame(my_data)
buys = data.loc[(data[6] == 1)]
sells = data.loc[(data[7]==-1)]

plt.figure()
plt.plot(data.index, data[3], label='Euro')
plt.plot(buys.index, data.iloc[buys.index][3], '^', markersize=10, color='g')
plt.plot(sells.index, data.iloc[sells.index][3], 'v', markersize=10, color='r')


my_data = TD_reverse_differential(my_data,4,5,8,9)


data = pd.DataFrame(my_data)
buys = data.loc[(data[8] == 1)]
sells = data.loc[(data[9]==-1)]
plt.figure()
plt.plot(data.index, data[3], label='Euro')
plt.plot(buys.index, data.iloc[buys.index][3], '^', markersize=10, color='g')
plt.plot(sells.index, data.iloc[sells.index][3], 'v', markersize=10, color='r')


my_data = TD_anti_differential(my_data,4,5,10,11)


data = pd.DataFrame(my_data)
buys = data.loc[(data[10] == 1)]
sells = data.loc[(data[11]==-1)]
plt.figure()
plt.plot(data.index, data[3], label='Euro')
plt.plot(buys.index, data.iloc[buys.index][3], '^', markersize=10, color='g')
plt.plot(sells.index, data.iloc[sells.index][3], 'v', markersize=10, color='r')


