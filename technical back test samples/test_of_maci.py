# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 13:42:46 2021

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


def jump(Data, jump):
    
    Data = Data[jump:, ]
    
    return Data
############
def ma(Data, lookback, what, where):
    
    for i in range(len(Data)):
            try:
                Data[i, where] = (Data[i - lookback + 1:i + 1, what].mean())
        
            except IndexError:
                pass
    return Data
############
def normalizer(Data, lookback, what, where):
        
    for i in range(len(Data)):
        
        try:
            Data[i, where] = (Data[i, what] - min(Data[i - lookback + 1:i + 1, what])) / (max(Data[i - lookback + 1:i + 1, what]) - min(Data[i - lookback + 1:i + 1, what]))
        
        except ValueError:
            pass
    
    Data[:, where] = Data[:, where] * 100    
    Data = jump(Data, lookback)
    return Data
############
def maci(Data, lookback, normalization_lookback, what, where):
    
    Data = ma(Data, lookback, what, where)
    
    Data[:, where + 1] = Data[:, what] - Data[:, where]
    
    Data = jump(Data, lookback)
    
    Data = normalizer(Data, normalization_lookback, 6, 7)
    
    Data = deleter(Data, where, 2)
    return Data

'''

Creating the Strategy
Now we can create a simple contrarian strategy based on the MACI(100, 10). This nomenclature means that it is a 100-period moving average with a normalization period of 10. The strategy is already obvious. We will use:
The 0 level of the MACI to trigger a long (buy) order and hold it until getting another signal or exiting following the Simplified ATR.
The 100 level of the MACI to trigger a short (sell) order and hold it until getting another signal or exiting following the Simplified ATR.


'''

