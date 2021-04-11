# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 12:34:59 2021

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
from various_ma_package import ema

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

def atr(Data, atr_lookback, high, low, close, whereTR, whereATR):
    
    # TR
    for i in range(len(Data)):
            
        Data[i, whereTR] = max(Data[i, high] - Data[i, low],
                    abs(Data[i, high] - Data[i - 1, close]),
                    abs(Data[i, low] - Data[i - 1, close]))
        
    Data[0, whereATR] = Data[0, high] - Data[0, low]
        
    Data[atr_lookback - 1, whereATR] = Data[0:atr_lookback - 1,  whereTR].mean()
    
    first_atr = Data[0, whereATR]
    
    for i in range(len(Data)):
        
        Data[0, whereATR] = first_atr 
        
        Data[i, whereATR] = (((Data[i - 1, whereATR] * (atr_lookback - 1)) + Data[i, whereTR]) / atr_lookback)
    
    return Data


def eATR(Data, atr_lookback, high, low, close, whereTR):
    for i in range(len(Data)):
        try:
            Data[i, whereTR] = max(Data[i,high] - Data[i, low], abs(Data[i, high] - Data[i-1, close]), abs(Data[i, low] - Data[i-1, close]))
        except ValueError:
            pass
        Data[0, whereTR] = 0 
        Data = ema(Data, 2, atr_lookback, whereTR, whereTR + 1, whereTR + 2)
    return Data

def supertrend(Data, multiplier, lookback):
    
    for i in range(len(Data)):
        
            # Average Price
            Data[i, 5] = (Data[i, 1] + Data[i, 2]) / 2
            # Basic Upper Band
            Data[i, 6] = Data[i, 5] + (multiplier * Data[i, 4])
            # Lower Upper Band
            Data[i, 7] = Data[i, 5] - (multiplier * Data[i, 4])
    
    # Final Upper Band
    for i in range(len(Data)):
        
        if i == 0:
            Data[i, 8] = 0
            
        else:  
            if (Data[i, 6] < Data[i - 1, 8]) or (Data[i - 1, 3] > Data[i - 1, 8]):
                Data[i, 8] = Data[i, 6]
            
            else:
                Data[i, 8] = Data[i - 1, 8]
    
    # Final Lower Band
    for i in range(len(Data)):
        
        if i == 0:
            Data[i, 9] = 0
            
        else:  
            if (Data[i, 7] > Data[i - 1, 9]) or (Data[i - 1, 3] < Data[i - 1, 9]):
                Data[i, 9] = Data[i, 7]
            
            else:
                Data[i, 9] = Data[i - 1, 9]
      
    # SuperTrend
    for i in range(len(Data)):
        
        if i == 0:
            Data[i, 10] = 0
        
        elif (Data[i - 1, 10] == Data[i - 1, 8]) and (Data[i, 3] <= Data[i, 8]):
            Data[i, 10] = Data[i, 8]
        
        elif (Data[i - 1, 10] == Data[i - 1, 8]) and (Data[i, 3] > Data[i, 8]):
            Data[i, 10] = Data[i, 9]
        
        elif (Data[i - 1, 10] == Data[i - 1, 9]) and (Data[i, 3] >= Data[i, 9]):
            Data[i, 10] = Data[i, 9]
        
        elif (Data[i - 1, 10] == Data[i - 1, 9]) and (Data[i, 3] < Data[i, 9]):
            Data[i, 10] = Data[i, 8]        
    
    return Data
