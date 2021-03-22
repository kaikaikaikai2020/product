# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 08:02:47 2021

@author: ASUS
"""
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



def stochastic(Data, lookback, onwhat, where):
        
    for i in range(len(Data)):
        
        try:
    
            Data[i, where] = (Data[i, onwhat] - min(Data[i -  lookback:i + 1, 2])) / (max(Data[i - lookback:i + 1, 1]) - min(Data[i - lookback:i + 1, 2]))
        except ValueError:
            pass
    
    Data[:, where] = Data[:, where] * 100
    return Data

def hidden_divergence(Data, lower_barrier, upper_barrier, width):
    
    for i in range(len(Data)):
        
        try:
            if Data[i, 4] < lower_barrier:
                
                for a in range(i + 1, i + width):
                    
                    
                    if Data[a, 4] > lower_barrier:
                        
                        for r in range(a + 1, a + width):
                            
                            if Data[r, 4] < lower_barrier and \
                            Data[r, 4] < Data[i, 4] and Data[r, 3] > Data[i, 3]:
                                
                                for s in range(r + 1, r + width):
                                    
                                   
                                    if Data[s, 4] > lower_barrier:
                                        Data[s, 6] = 1
                                        break
                                    
                                    else:
                                        break
                            else:
                                break
                        else:
                            break
                    else:
                        break
                    
        except IndexError:
            pass
        
    for i in range(len(Data)):
        
        try:
            if Data[i, 4] > upper_barrier:
                
                for a in range(i + 1, i + width):
                    
                    
                    if Data[a, 4] < upper_barrier:
                        for r in range(a + 1, a + width):
                            if Data[r, 4] > upper_barrier and \
                            Data[r, 4] > Data[i, 4] and Data[r, 3] < Data[i, 3]:
                                for s in range(r + 1, r + width):
                                    
                                    
                                    if Data[s, 4] < upper_barrier:
                                        Data[s, 7] = -1
                                        break
                                    else:
                                        break
                            else:
                                break
                        else:
                            break
                    else:
                        break
        except IndexError:
            pass 
        
    return Data

def adder(Data, times):
    for i in range(1, times + 1):
        z = np.zeros((len(Data), 1), dtype = float)
        Data = np.append(Data, z, axis = 1)
    return Data



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
    my_data = stochastic(my_data, lookback,3,4)
    my_data = hidden_divergence(my_data, low_barrier, up_barrier,width)
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
