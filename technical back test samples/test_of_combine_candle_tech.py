# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 11:40:35 2021

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

my_data = np.round(my_data,4)

def doji_scan(Data, buy, sell):
    for i in range(len(Data)):
        
        if Data[i - 1, 3] < Data[i - 1, 0 ] and Data[i, 3] == Data[i, 0]:
            Data[i, buy] = 1
        
        if Data[i - 1, 3] > Data[i - 1, 0 ] and Data[i, 3] == Data[i, 0]:
            Data[i, sell] = -1
            
    return Data


# The variable Data refers to the OHLC Historical Data array.
# The variable buy refers to where the long order is triggered
# The variable sell refers to where the short order is triggered


# need to calculate RSI first


def complex_strategy_1(Data, buy, sell):
 for i in range(len(Data)):
        
  if Data[i, rsi] < lower and Data[i, 3] == Data[i, 0]:
   Data[i, buy] = 1
        
  if Data[i, rsi] > upper and Data[i, 3] == Data[i, 0]:
   Data[i, sell] = -1
# The variable Data refers to the OHLC Historical Data array.
# The variable indicator refers to RSI's location.
# The variable buy refers to where the long order is triggered.
# The variable sell refers to where the short order is triggered.

def harami_scan(Data, buy, sell):
    for i in range(len(Data)):
        if Data[i - 1, 3] < Data[i - 1, 0] and Data[i, 3] > Data[i, 0] and Data[i - 1, 3] < Data[i, 2] and  Data[i - 1, 0] > Data[i, 1]:
            Data[i, buy] = 1
        
        if Data[i - 1, 3] > Data[i - 1, 0] and Data[i, 3] < Data[i, 0] and Data[i - 1, 3] > Data[i, 1] and  Data[i - 1, 0] < Data[i, 2]:
            Data[i, sell] = -1
     return Data

# The variable Data refers to the OHLC Historical Data array.
# The variable buy refers to where the long order is triggered
# The variable sell refers to where the short order is triggered
# Defining the barriers
lower = 45
upper = 55
def complex_strategy_2(Data, indicator, buy, sell):
    for i in range(len(Data)):
        if Data[i, indicator] < lower and Data[i - 1, 3] < Data[i - 1, 0] and  Data[i, 3] > Data[i, 0] and Data[i - 1, 3] < Data[i, 2] and  Data[i - 1, 0] > Data[i, 1]:
            Data[i, buy] = 1
        
        if Data[i, indicator] > upper and Data[i - 1, 3] > Data[i - 1, 0] and Data[i, 3] < Data[i, 0] and Data[i - 1, 3] > Data[i, 1] and Data[i - 1, 0] < Data[i, 2]:
            Data[i, sell] = -1
            
    return Data
# The variable Data refers to the OHLC Historical Data array.
# The variable indicator refers to RSI's location.
# The variable buy refers to where the long order is triggered.
# The variable sell refers to where the short order is triggered.

def piercing_scan(Data, buy, sell):
    for i in range(len(Data)):
        
        if  Data[i - 1, 3] < Data[i - 1, 0] and \
            Data[i, 3] > Data[i, 0] and \
            Data[i, 0] < Data[i - 1, 3] and \
            Data[i, 3] > Data[i - 1, 3] and \
            Data[i, 3] < Data[i - 1, 0]:
            Data[i, buy] = 1
            
    return Data

def cloud_scan(Data, buy, sell):
    for i in range(len(Data)):
        if  Data[i - 1, 3] > Data[i - 1, 0] and \
            Data[i, 3] < Data[i, 0] and \
            Data[i, 0] > Data[i - 1, 3] and \
            Data[i, 3] < Data[i - 1, 3] and \
            Data[i, 3] > Data[i - 1, 0]:
            Data[i, sell] = -1
    return Data

        
# Normally, we should put them together inside one function.
# The variable Data refers to the OHLC Historical Data array.
# The variable buy refers to where the long order is triggered
# The variable sell refers to where the short order is triggered
def complex_strategy_3(Data, buy, sell):
    for i in range(len(Data)):
        
        if Data[i, indicator] < lower and \
            Data[i - 1, 3] < Data[i - 1, 0] and \
            Data[i, 3] > Data[i, 0] and \
            Data[i, 0] < Data[i - 1, 3] and \
            Data[i, 3] > Data[i - 1, 3] and \
            Data[i, 3] < Data[i - 1, 0]:
            Data[i, buy] = 1
        
        if Data[i, indicator] > upper and \
            Data[i - 1, 3] > Data[i - 1, 0] and \
            Data[i, 3] < Data[i, 0] and \
            Data[i, 0] > Data[i - 1, 3] and \
            Data[i, 3] < Data[i - 1, 3] and \
            Data[i, 3] > Data[i - 1, 0]:
            Data[i, sell] = -1
            
    return Data
# The variable Data refers to the OHLC Historical Data array.
# The variable indicator refers to RSI's location.
# The variable buy refers to where the long order is triggered.
# The variable sell refers to where the short order is triggered.