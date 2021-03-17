# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 13:16:17 2021

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
lookback = 5

#fld_list = ['open', 'close', 'high', 'low']
#prev = (today - BDay(5))
#my_data = bbg.bdib(ticker, fld_list, datetime(2021,1,1,1,30), datetime(2021,2,1,21,30), eventType='TRADE', interval = 30)
my_data=bbg.bdh(ticker, fld_list, '20131201', end_date=date.today().strftime('%Y%m%d'), periodselection = 'DAILY')
#my_data = bbg.bdib(ticker, fld_list, prev, today, eventType='TRADE', interval = 60)
print(my_data)

def momentum_indicator(Data, what, where, lookback):
    
 for i in range(len(Data)):
     
     Data[i, where] = Data[i, what] / Data[i - lookback, what] * 100
    
 return Data
# The Data variable is the OHLC array
# The what variable is the closing price column
# The where variable is where to put the indicator
# The lookback variable is the subtraction range


def astral(Data, completion, step, step_two, what, high, low, where_long, where_short):
   
    # Timing buy signal
    counter = -1
    for i in range(len(Data)):    
        if Data[i, what] < Data[i - step, what] and Data[i, low] < Data[i - step_two, low]:
            Data[i, where_long] = counter
            counter += -1       
            if counter == -completion - 1:
                counter = 0
            else:
                continue        
        elif Data[i, what] >= Data[i - step, what]:
            counter = -1 
            Data[i, where_long] = 0 
        
    # Timing sell signal       
    counter = 1 
    
    for i in range(len(Data)):
        if Data[i, what] > Data[i - step, what] and Data[i, high] > Data[i - step_two, high]: 
            Data[i, where_short] = counter 
            counter += 1        
            if counter == completion + 1: 
                counter = 0            
            else:
                continue        
        elif Data[i, what] <= Data[i - step, what]: 
            counter = 1 
            Data[i, where_short] = 0 
      
    return Data

# The completion variable refers to Astral's final count
# The step variable refers to Astral's first lookback
# The step_two variable refers to Astral's second lookback
# The what variable refers to the closing price
# The high variable refers to the high price
# The low variable refers to the low price
# The where_long variable refers to where to put the buy trigger
# The where_short variable refers to where to put the sell trigger


def adder(Data, times):
    for i in range(1, times + 1):
        z = np.zeros((len(Data), 1), dtype = float)
        Data = np.append(Data, z, axis = 1)
    return Data
# Using the function to add 20 columns
my_data = adder(my_data, 20)

#my_data = momentum_indicator(my_data, lookback, 4, 5)
my_data = momentum_indicator(my_data, what=3, where=4, lookback=5)
my_data = astral(my_data, 8, 3, 5, 3, 1, 2, 6, 7)

def simplified_astral(Data, completion, step, what, where_long, where_short):
        
    # Timing buy signal
    counter = -1
    for i in range(len(Data)):    
        if Data[i, what] < Data[i - step, what]:
            Data[i, where_long] = counter
            counter += -1       
            if counter == -completion - 1:
                counter = 0
            else:
                continue        
        elif Data[i, what] >= Data[i - step, what]:
            counter = -1 
            Data[i, where_long] = 0 
    
    if Data[8, where_long] == -completion:
        Data = Data[9:,]
    elif Data[7, where_long] == -completion + 1:
        Data = Data[8:,]
    elif Data[6, where_long] == -completion + 2:
        Data = Data[7:,]
    elif Data[5, where_long] == -completion + 3:
        Data = Data[6:,]
    elif Data[4, where_long] == -completion + 4:
        Data = Data[5:,]
        
    # Timing sell signal       
    counter = 1 
    
    for i in range(len(Data)):
        if Data[i, what] > Data[i - step, what]: 
            Data[i, where_short] = counter 
            counter += 1        
            if counter == completion + 1: 
                counter = 0            
            else:
                continue        
        elif Data[i, what] <= Data[i - step, what]: 
            counter = 1 
            Data[i, where_short] = 0 
    
    if Data[8, where_short] == completion:
        Data = Data[9:,]
    elif Data[7, where_short] == completion - 1:
        Data = Data[8:,]
    elif Data[6, where_short] == completion - 2:
        Data = Data[7:,]
    elif Data[5, where_short] == completion - 3:
        Data = Data[6:,]
    elif Data[4, where_short] == completion - 4:
        Data = Data[5:,] 
        
    return Data
# Using the function on the Momentum Indicator
lower_barrier =-8
upper_barrier = 8
completion = 8
step =3 

my_data = simplified_astral(my_data, completion, step, 5, 8, 9)
my_data[:, 10] = my_data[:, 6] + my_data[:, 8]
my_data[:, 11] = my_data[:, 7] + my_data[:, 9]


def signal(Data, astral_buy, astral_sell, buy, sell, lower_barrier, upper_barrier):
    
    for i in range(len(Data)):
            
        if Data[i, astral_buy] <= lower_barrier and Data[i - 1, astral_buy] > lower_barrier:
            Data[i, buy] = 1
            
        if Data[i, astral_sell] >= upper_barrier and Data[i - 1, astral_sell] < upper_barrier:
            Data[i, sell] = -1
            
    return Data
            
#print(my_data)        
my_data = signal(my_data,6,7,12,13, lower_barrier, upper_barrier)

data = pd.DataFrame(my_data)
buys = data.loc[(data[12] == 1)]
sells = data.loc[(data[13]==-1)]

plt.plot(data.index, data[3], label='Euro')
plt.plot(buys.index, data.iloc[buys.index][3], '^', markersize=10, color='g')
plt.plot(sells.index, data.iloc[sells.index][3], 'v', markersize=10, color='r')