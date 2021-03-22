
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


for ticker in instrument_list:
    #fld_list = ['open', 'close', 'high', 'low']
    #prev = (today - BDay(5))
    #my_data = bbg.bdib(ticker, fld_list, datetime(2021,1,1,1,30), datetime(2021,2,1,21,30), eventType='TRADE', interval = 30)
    my_data=bbg.bdh(ticker, fld_list, '20131201', end_date=date.today().strftime('%Y%m%d'), periodselection = 'DAILY')
    #my_data = bbg.bdib(ticker, fld_list, prev, today, eventType='TRADE', interval = 60)
    print(my_data)
    lookback = 14
    
    Data = my_data
    delta = Data.iloc[:, 3].diff()
    delta = delta[1:]
    # Make the positive gains (up) and negative gains (down) Series
    up, down = delta.copy(), delta.copy()
    up[up < 0] = 0
    down[down > 0] = 0
    roll_up = up.ewm(span=lookback).mean()
    roll_down = (down.abs().ewm(span=lookback)).mean()
    
    roll_up = roll_up[lookback:]
    roll_down = roll_down[lookback:]
    Data = Data.iloc[lookback + 1:,].values
    # Calculate the SMA
    # Calculate the RSI based on SMA
    RS = roll_up / roll_down
    RSI = (100.0 - (100.0 / (1.0 + RS)))
    RSI = np.array(RSI)
    RSI = np.reshape(RSI, (-1, 1))
    Data = np.concatenate((Data, RSI), axis = 1)
    
    def adder(Data, times):
        for i in range(1, times + 1):
            z = np.zeros((len(Data), 1), dtype = float)
            Data = np.append(Data, z, axis = 1)
        return Data
    # Using the function to add 20 columns
    Data = adder(Data, 20)
    
    
    lower_barrier = 30
    upper_barrier = 70
    width = 10
    
    
    # Bullish Divergence
    for i in range(len(Data)):
        try:
            if Data[i, 4] < lower_barrier:
                
               for a in range(i + 1, i + width):
                   if Data[a, 4] > lower_barrier:
                      
                      for r in range(a + 1, a + width):
                         if Data[r, 4] < lower_barrier and Data[r, 4] > Data[i, 4] and Data[r, 3] < Data[i, 3]:
                               
                              for s in range(r + 1, r + width): 
                                  if Data[s, 4] > lower_barrier:
                                      Data[s + 1, 5] = 1
                                      break
                                    
                                  else:
                                      continue
                         else:
                              continue
                   else:
                           continue
            else:
                continue
                    
        except IndexError:
             pass
    
    # Bearish Divergence
    for i in range(len(Data)):
       
        try:
            if Data[i, 4] > upper_barrier:
                
                for a in range(i + 1, i + width): 
                    if Data[a, 4] < upper_barrier:
                        for r in range(a + 1, a + width):
                            if Data[r, 4] > upper_barrier and \
                            Data[r, 4] < Data[i, 4] and Data[r, 3] >Data[i, 3]:
                                for s in range(r + 1, r + width):
                                    if Data[s, 4] < upper_barrier:
                                        Data[s + 1, 6] = -1
                                        break
                                    else:
                                        continue
                            else:
                                continue
                        else:
                            continue
                    else:
                        continue
        except IndexError:
            pass
    
    
    data = pd.DataFrame(Data)
    buys = data.loc[(data[5] == 1)]
    sells = data.loc[(data[6]==-1)]
    plt.figure()
    plt.plot(data.index, data[3], label='Euro')
    plt.plot(buys.index, data.iloc[buys.index][3], '^', markersize=10, color='g')
    plt.plot(sells.index, data.iloc[sells.index][3], 'v', markersize=10, color='r')
    
    my_data = Data
    my_data[:,8] = my_data[:,5]+my_data[:,6]
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
    
    
