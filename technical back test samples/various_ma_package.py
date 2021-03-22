# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 12:32:46 2021

@author: ASUS

"""
import numpy as np
import pandas as pd

def ma(Data, period, onwhat, where):
    
    for i in range(len(Data)):
            try:
                Data[i, where] = (Data[i - period:i + 1, onwhat].mean())
        
            except IndexError:
                pass
    return Data

def ema(Data, alpha, lookback, what, where):
    
    # alpha is the smoothing factor
    # window is the lookback period
    # what is the column that needs to have its average calculated
    # where is where to put the exponential moving average
    
    alpha = alpha / (lookback + 1.0)
    beta  = 1 - alpha
    
    # First value is a simple SMA
    Data = ma(Data, lookback, what, where)
    
    # Calculating first EMA
    Data[lookback + 1, where] = (Data[lookback + 1, what] * alpha) + (Data[lookback, where] * beta)
    # Calculating the rest of EMA
    for i in range(lookback + 2, len(Data)):
            try:
                Data[i, where] = (Data[i, what] * alpha) + (Data[i - 1, where] * beta)
        
            except IndexError:
                pass
    return Data

def atr(Data, lookback, high, low, close, where):
   # TR
   for i in range(len(Data)):
       try:
           Data[i, where] = max(Data[i, high] - Data[i, low],abs(Data[i, high] - Data[i - 1, close]), abs(Data[i, low] - Data[i - 1, close]))
       except ValueError:
            pass
   Data[0, where] = 0    
    # ATR
   Data = ma(Data, lookback, where, where + 1)
   return Data

def keltner_channel(Data, ma_lookback, atr_lookback, multiplier, what, where):
    Data = ema(Data, 2, ma_lookback, what, where)
    Data = atr(Data, atr_lookback, 2, 1, 3, where + 1)
    Data[:, where + 2] = Data[:, where] + (Data[:, where + 1] * multiplier)
    Data[:, where + 3] = Data[:, where] - (Data[:, where + 1] * multiplier)
    return Data

def ema_old(Data, alpha, window, what, whereSMA, whereEMA):
    
    # alpha is the smoothing factor
    # window is the lookback period
    # what is the column that needs to have its average calculated
    # where is where to put the exponential moving average
    
    alpha = alpha / (window + 1.0)
    beta  = 1 - alpha
    
    # First value is a simple SMA
    Data[window - 1, whereSMA] = np.mean(Data[:window - 1, what])
    
    # Calculating first EMA
    Data[window, whereEMA] = (Data[window, what] * alpha) + (Data[window - 1, whereSMA] * beta)
# Calculating the rest of EMA
    for i in range(window + 1, len(Data)):
            try:
                Data[i, whereEMA] = (Data[i, what] * alpha) + (Data[i - 1, whereEMA] * beta)
        
            except IndexError:
                pass
    return Data




def lwma(Data, period):
    weighted = []
    for i in range(len(Data)):
            try:
                total = np.arange(1, period + 1, 1) # weight matrix
                
                matrix = Data[i - period + 1: i + 1, 3:4]
                matrix = np.ndarray.flatten(matrix)
                matrix = total * matrix # multiplication
                wma = (matrix.sum()) / (total.sum()) # WMA
                weighted = np.append(weighted, wma) # add to array
            except ValueError:
                pass
    return weighted


def kama(Data, what, where, change):
    
    # Change from previous period
    for i in range(len(Data)):
        Data[i, where] = abs(Data[i, what] - Data[i - 1, what])
    
    Data[0, where] = 0
    
    # Sum of changes
    for i in range(len(Data)):
        Data[i, where + 1] = (Data[i - change + 1:i + 1, where].sum())   
        
    # Volatility    
    for i in range(len(Data)):
        Data[i, where + 2] = abs(Data[i, 3] - Data[i - 10, 3])
        
    Data = Data[11:, ]
    
    # Efficiency Ratio
    Data[:, where + 3] = Data[:, where + 2] / Data[:, where + 1]
    
    for i in range(len(Data)):
        Data[i, where + 4] = np.square(Data[i, where + 3] * 0.6666666666666666667)
        
    for i in range(len(Data)):
        Data[i, where + 5] = Data[i - 1, where + 5] + (Data[i, where + 4] * (Data[i, 3] - Data[i - 1, where + 5]))
        Data[11, where + 5] = 0

def adder(Data, times):
    for i in range(1, times + 1):
        z = np.zeros((len(Data), 1), dtype = float)
        Data = np.append(Data, z, axis = 1)
    return Data

def rsi(Data, rsi_lookback, what1, what2):
    
    rsi_lookback = (rsi_lookback * 2) - 1 # From exponential to smoothed
          
    # Get the difference in price from previous step
    delta = []
   
    for i in range(len(Data)):
        try:
            diff = Data[i, what1] - Data[i - 1, what1] 
            delta = np.append(delta, diff)                  
        except IndexError:
            pass
        
    delta = np.insert(delta, 0, 0, axis = 0)               
    delta = delta[1:] 
    
    # Make the positive gains (up) and negative gains (down) Series
    up, down = delta.copy(), delta.copy()
    up[up < 0] = 0
    down[down > 0] = 0
    
    up = np.array(up)
    down = np.array(down)
    
    roll_up = up
    roll_down = down
    
    roll_up = np.reshape(roll_up, (-1, 1))
    roll_down = np.reshape(roll_down, (-1, 1))
    
    roll_up = adder(roll_up, 3)
    roll_down = adder(roll_down, 3)
    
    roll_up = ema(roll_up, 2, rsi_lookback, what2, 1)
    roll_down = ema(abs(roll_down), 2, rsi_lookback, what2, 1)
    
    # Calculate the SMA
    roll_up = roll_up[rsi_lookback:, 1:2]
    roll_down = roll_down[rsi_lookback:, 1:2]
    Data = Data[rsi_lookback + 1:,]
    
    # Calculate the RSI based on SMA
    RS = roll_up / roll_down
    RSI = (100.0 - (100.0 / (1.0 + RS)))
    RSI = np.array(RSI)
    RSI = np.reshape(RSI, (-1, 1))
    RSI = RSI[1:,]
    
    Data = np.concatenate((Data, RSI), axis = 1)
    return Data


def momentum_indicator(Data, what, where, lookback):
    for i in range(len(Data)):
     Data[i, where] = Data[i, what] / Data[i - lookback, what] * 100
    
    return Data
# The Data variable is the OHLC array
# The what variable is the closing price column
# The where variable is where to put the indicator
# The lookback variable is the subtraction range which is 5

from scipy.stats import pearsonr

def jump(Data, jump):
    
    Data = Data[jump:, ]
    
    return Data



def rolling_correlation(Data, first_data, second_data, lookback, where):
    
    for i in range(len(Data)):
        
        try:
            Data[i, where] = pearsonr(Data[i - lookback + 1:i + 1, first_data], Data[i - lookback + 1:i + 1, second_data])[0]
            
            
        except ValueError:
            pass
    
    Data = jump(Data, lookback) 
    
    return Data

def normalizer(Data, lookback, what, where):
        
    for i in range(len(Data)):
        
        try:
            Data[i, where] = (Data[i, what] - min(Data[i - lookback + 1:i + 1, what])) / (max(Data[i - lookback + 1:i + 1, what]) - min(Data[i - lookback + 1:i + 1, what]))
        
        except ValueError:
            pass
    
    Data[:, where] = Data[:, where] * 100    
    Data = jump(Data, lookback)
    return Data


def stochastics(Data, lookback, what, high, low, where):
        
    for i in range(len(Data)):
        
        try:
            Data[i, where] = (Data[i, what] - min(Data[i - lookback + 1:i + 1, low])) / (max(Data[i - lookback + 1:i + 1, high]) - min(Data[i - lookback + 1:i + 1, low]))
        
        except ValueError:
            pass
    
    Data[:, where] = Data[:, where] * 100
    return Data
# The Data variable refers to the OHLC array
# The lookback variable refers to the period (5, 14, 21, etc.)
# The what variable refers to the closing price
# The high variable refers to the high price
# The low variable refers to the low price
# The where variable refers to where to put the Oscillator



def volatility(Data, lookback, what, where):
    
    for i in range(len(Data)):
            try:
                Data[i, where] = (Data[i - lookback + 1:i + 1, what].std())
        
            except IndexError:
                pass
    return Data


'''

We can use more complex stuff to account for the lookbacks but we can try something more creative and use insights from the Fibonacci sequence. Here is how:
Whenever the normalized values of volatility are between 0 and 10, the lookback period of the Stochastics is 144.
Whenever the normalized values of volatility are between 10 and 20, the lookback period of the Stochastics is 89.
Whenever the normalized values of volatility are between 20 and 30, the lookback period of the Stochastics is 55.
Whenever the normalized values of volatility are between 30 and 40, the lookback period of the Stochastics is 34.
Whenever the normalized values of volatility are between 40 and 50, the lookback period of the Stochastics is 21.
Whenever the normalized values of volatility are between 50 and 60, the lookback period of the Stochastics is 13.
Whenever the normalized values of volatility are between 60 and 70, the lookback period of the Stochastics is 8.
Whenever the normalized values of volatility are between 70 and 80, the lookback period of the Stochastics is 5.
Whenever the normalized values of volatility are between 80 and 90, the lookback period of the Stochastics is 3.
Whenever the normalized values of volatility are between 90 and 100, the lookback period of the Stochastics is 2.

'''

def fib_stoch(Data, volatility_lookback, what, where):
        
    Data = volatility(Data, volatility_lookback, what, where)
    Data = normalizer(Data, volatility_lookback, where, where + 1)
    for i in range(len(Data)):
        Data[i, where + 1] = round(Data[i, where + 1], 0) 
    
    for i in range(len(Data)):
        
        if Data[i, where + 1] >= 0 and Data[i, where + 1] <= 10 :
            Data[i, where + 1] = 144
        if Data[i, where + 1] > 10 and Data[i, where + 1] <= 20 :
            Data[i, where + 1] = 89            
        if Data[i, where + 1] > 20 and Data[i, where + 1] <= 30 :
            Data[i, where + 1] = 55
        if Data[i, where + 1] > 30 and Data[i, where + 1] <= 40 :
            Data[i, where + 1] = 34
        if Data[i, where + 1] > 40 and Data[i, where + 1] <= 50 :
            Data[i, where + 1] = 21
        if Data[i, where + 1] > 50 and Data[i, where + 1] <= 60 :
            Data[i, where + 1] = 13
        if Data[i, where + 1] > 60 and Data[i, where + 1] <= 70 :
            Data[i, where + 1] = 8
        if Data[i, where + 1] > 70 and Data[i, where + 1] <= 80 :
            Data[i, where + 1] = 5
        if Data[i, where + 1] > 80 and Data[i, where + 1] <= 90 :
            Data[i, where + 1] = 3
        if Data[i, where + 1] > 90 and Data[i, where + 1] <= 100 :
            Data[i, where + 1] = 2            
    
    Data = jump(Data, volatility_lookback)
        
    for i in range(len(Data)):
        
        try:
            lookback = int(Data[i, where + 1])
            Data[i, where + 2] = (Data[i, what] - min(Data[i - lookback + 1:i + 1, 2])) / (max(Data[i - lookback + 1:i + 1, 1]) -    min(Data[i - lookback + 1:i + 1, 2]))
        
        except ValueError:
            pass
    return Data


#def ema(Data, alpha, window, what, whereSMA, whereEMA):
def fractal_indicator(Data, high, low, ema_lookback, min_max_lookback, where):
    
    Data = ema(Data, 2, ema_lookback, high, 20, where)
    Data = ema(Data, 2, ema_lookback, low, 20, where + 1)
    
    Data = volatility(Data, ema_lookback, high, where + 2)
    Data = volatility(Data, ema_lookback, low, where + 3)
    
    Data[:, where + 4] = Data[:, high] - Data[:, where]
    Data[:, where + 5] = Data[:, low]  - Data[:, where + 1]
    for i in range(len(Data)):
        try:
            Data[i, where + 6] = max(Data[i - min_max_lookback + 1:i + 1, where + 4])
        
        except ValueError:
            pass              
        
    for i in range(len(Data)):
        try:
            Data[i, where + 7] = min(Data[i - min_max_lookback + 1:i + 1, where + 5])
        
        except ValueError:
            pass  
         
    Data[:, where + 8] =  (Data[:, where +  2] + Data[:, where +  3]) / 2
    Data[:, where + 9] = (Data[:, where + 6] - Data[:, where + 7]) / Data[:, where + 8]
    
    return Data


def BollingerBands(Data, boll_lookback, standard_distance, what, where):
       
    # Calculating mean
    ma(Data, boll_lookback, what, where)
    
    # Calculating volatility
    volatility(Data, boll_lookback, what, where + 1)
    
    # Calculating Bands
    Data[:, where + 2] = Data[:, where] + (standard_distance * Data[:, where + 1])
    Data[:, where + 3] = Data[:, where] - (standard_distance * Data[:, where + 1])
        
    return Data

def adx(Data, high, low, close, lookback, where):
    # DM+
    for i in range(len(Data)):
        if (Data[i, high] - Data[i - 1, high]) > (Data[i - 1, low] - Data[i, low]):
            Data[i, where] = Data[i, high] - Data[i - 1, high]
            
        else:
            Data[i, where] = 0
        
    # DM-
    for i in range(len(Data)):
        if (Data[i, high] - Data[i - 1, high]) < (Data[i - 1, low] - Data[i, low]):
            Data[i, where + 1] = Data[i - 1, low] - Data[i, low]
            
        else:
            Data[i, where + 1] = 0
        
    # Smoothing DI+   
    Data = ema(Data, 2, (lookback * 2 - 1), where, where + 2)
    
    # Smoothing DI-
    Data = ema(Data, 2, (lookback * 2 - 1), where + 1, where + 3)
    
    # Smoothing ATR
    Data = atr(Data, (lookback * 2 - 1), high, low, close, where + 4)
    
    # DI+
    Data[:, where + 5] = Data[:, where + 2] / Data[:, where + 4]
    
    # DI-
    Data[:, where + 6] = Data[:, where + 3] / Data[:, where + 4]
    
    # ADX
    for i in range(len(Data)):
        Data[i, where + 7] = abs(Data[i, where + 5] - Data[i, where + 6]) / abs(Data[i, where + 5] + Data[i, where + 6]) * 100
    
    Data = ema(Data, 2, (lookback * 2 - 1), where + 7, where + 8)  
    
    return Data

def mawi(Data, short_ma, long_ma, normalization_lookback, what, where):
    
    Data = ma(Data, short_ma, what, where)
    Data = ma(Data, long_ma, what, where + 1)
    # MAWI line (Width) 
    Data[:, where + 2] = Data[:, where] - Data[:, where + 1]
    # MAWI normalized
    Data = normalizer(Data, normalization_lookback, where + 2, where + 3) 
    Data[:, where + 3] =  Data[:, where + 3]
    
    return Data

def rri(Data, lookback, where):
    
    for i in range(len(Data)):
        
        Data[i, where] = (Data[i, 3] - Data[i - lookback, 0]) / (Data[i - lookback, 1] - Data[i - lookback, 2])
        if Data[i - lookback, 1] == Data[i - lookback, 2]:
            Data[i, where] = 0
    
    return Data

def deleter(Data, index, times):
    
    for i in range(1, times + 1):
        Data = np.delete(Data, index, axis = 1)
    return Data

def cci(Data, lookback, where, constant):
    
    # Calculating Typical Price
    Data[:, where] = (Data[:, 1] + Data[:, 2] + Data[:, 3]) / 3
    
    # Calculating the Absolute Mean Deviation
    specimen = Data[:, where]
    MAD_Data = pd.Series(specimen)
    
    for i in range(len(Data)):
            
            Data[i, where + 1] = MAD_Data[i - lookback:i].mad()
    
    # Calculating Mean of Typical Price 
    Data = ma(Data, lookback, where, where + 2)
    
    # CCI
    for i in range(len(Data)):
      Data[i, where + 3] = (Data[i, where] - Data[i, where + 2]) / (constant * Data[i, where + 1]) 
    
    Data = jump(Data, lookback)   
    #Data = deleter(Data, where, 3)
    
    return Data

def countdown_indicator(Data, lookback, ma_lookback, opening, high, low, close, where):
    
    # Adding columns
    Data = adder(Data, 20)
    
    # Calculating Upside Pressure
    for i in range(len(Data)):
        
        if Data[i, close] > Data[i, opening]:
            
            Data[i, where] = 1
        if Data[i, high] > Data[i - 1, high]:
            
            Data[i, where + 1] = 1
            
    Data[:, where + 2] = Data[:, where] + Data[:, where + 1]    
    Data = deleter(Data, where, 2)    
            
    # Calculating Downside Pressure
    for i in range(len(Data)):
        
        if Data[i, close] < Data[i, opening]:
            
            Data[i, where + 1] = 1
        if Data[i, low] < Data[i - 1, low]:
            
            Data[i, where + 2] = 1
            
    Data[:, where + 3] = Data[:, where + 1] + Data[:, where + 2]    
    Data = deleter(Data, where + 1, 2)     
    
    # Calculate Cumulative Upside Pressure
    for i in range(len(Data)):
        
        Data[i, where + 2] = Data[i - lookback + 1:i + 1, where].sum()
    
    # Calculate Cumulative Downside Pressure
    for i in range(len(Data)):
        
        Data[i, where + 3] = Data[i - lookback + 1:i + 1, where + 1].sum()
       
    # Calculate the Countdown Indicator
    Data[:, where + 4] = Data[:, where + 2] - Data[:, where + 3]    
    Data = ema(Data, 2, ma_lookback, where + 4, where + 5)
    
    Data = deleter(Data, where, 5)
    Data = jump(Data, lookback)
       
    return Data