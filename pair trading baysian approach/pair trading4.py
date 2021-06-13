# -*- coding: utf-8 -*-
"""
Created on Sat Jun 12 21:31:24 2021

@author: ASUS
"""
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 12 20:06:43 2021

@author: ASUS
"""
#import pybbg 
import pandas as pd
from datetime import datetime, date, time
from pandas.tseries.offsets import BDay
import numpy as np
from scipy.signal import argrelextrema

today = datetime.today()

#prev = (today - BDay(10))
prev = (today - BDay(5))


result = pd.DataFrame()
from tqdm import tqdm
import matplotlib.pyplot as plt
#col_list = ["Date", "MSFT US EQUITY","AAPL US EQUITY"]

#stock1_name = '8801 JT EQUITY'
#stock2_name = '8802 JT EQUITY'

#stock1_name = 'HSI INDEX'
#stock2_name = 'HSCEI INDEX'



#stock1_name = '601318 CH EQUITY'
#stock2_name = '2318 CH EQUITY'

#stock1_name = 'HSCEI INDEX'
#stock2_name = 'XIN9I INDEX'

stock1_name = 'MSFT US EQUITY'
stock2_name = 'AAPL US EQUITY'

result = pd.read_csv('temp_data10.csv',dtype={stock1_name:float,stock2_name:float})    
#result.rename(columns=col_list, inplace=True)
result = result.set_index('Date')
print(result)

print('########## End process TW Stocks ########')

#result.plot()
#plt.show()
#data1 = result['MSFT US EQUITY'].to_numpy()
#data2 = result['AAPL US EQUITY'].to_numpy()


#data = orgin_data.diff().cumsum()

#print(result[stock1_name])
#orgin_data = result[[stock1_name,stock2_name]]
#result[stock1_name] = pd.to_numeric(result[stock1_name])
#result[stock2_name]= pd.to_numeric(result[stock2_name])
orgin_data = result
print(orgin_data.info())
#orgin_data[stock1_name]=orgin_data[stock1_name].astype(float)
#orgin_data[stock2_name]=orgin_data[stock2_name].astype(float)

#print(orgin_data)
print('----------------------------------------')
data = orgin_data.diff().cumsum()
#data = (orgin_data.pct_change()+1).cumprod()
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
#scaler = StandardScaler()
#scaler.fit(data)
#data = scaler.transform(data)
#print(data)
#data1 = data[0]
#data2 = data[1]
#print(data1)
#print(data2)


x = data.values #returns a numpy array
#min_max_scaler = MinMaxScaler()
sc = StandardScaler()
x_scaled = sc.fit_transform(x)
data = pd.DataFrame(x_scaled,columns=[stock1_name,stock2_name])
print(data)
data1 = data[stock1_name].ffill().fillna(0).values
data2 = data[stock2_name].ffill().fillna(0).values

#import statsmodels.api as sm
#model = sm.OLS(data2,data1)
#results = model.fit()
#print(results.params)
#mu = results.params
#mu = np.repeat(results.params, data1.shape[0])
#print(mu)

import requests

import pymc3 as pm
import theano as th
import seaborn as sns

import matplotlib.cm as cm
from matplotlib.collections import LineCollection

from sklearn.preprocessing import MinMaxScaler
from sklearn import cluster, covariance, manifold
from theano import shared

with pm.Model() as model:
    
    # inject external stock data
    stock1 = th.shared(data1)
    stock2 = th.shared(data2)
    
    # define our cointegration variables
    beta_sigma = pm.Exponential('beta_sigma', 50)
    beta = pm.GaussianRandomWalk('beta',sd=beta_sigma,
                                 shape=data1.shape[0])
    
    # with our assumptions, cointegration can be reframed as a regression problem
    stock2_regression = beta * stock1
    print(stock2_regression)

    # Assume prices are Normally distributed, the mean comes from the regression.
    sd = pm.HalfNormal('sd', sd=.1)
    likelihood = pm.Normal('y',
                           mu=stock2_regression,
                           sd=sd,
                           observed=stock2)
print(data1)
print(data2)

with model:
    stock1.set_value(data1)
    stock2.set_value(data2)
    #trace = pm.sample(2000,init="adapt_diag",tune=1000,cores=4)
    trace = pm.sample(2000,tune=1000,cores=4)

rolling_beta = trace[beta].T.mean(axis=1)

plt.figure(figsize = (18,8))
ax = plt.gca()
plt.title("Beta Distribution over Time")
pd.Series(rolling_beta).plot(ax=ax,color='r',zorder=1e6,linewidth=2)
for orbit in trace[beta][:500]:
    pd.Series(orbit).plot(ax=ax,color=sns.color_palette()[0],alpha=0.05)
plt.legend(['Beta Mean','Beta Orbit'])
#plt.savefig("beta distrib.png")

plt.show()
plt.savefig('1.png')

def getStrategyPortfolioWeights(rolling_beta,stock_name1,stock_name2,data,smoothing_window=15):

    data1 = data[stock_name1].ffill().fillna(0).values
    data2 = data[stock_name2].ffill().fillna(0).values

    # initial signal rebalance
    fixed_beta = rolling_beta[smoothing_window]
    signal = fixed_beta*data1 - data2
    smoothed_signal = pd.Series(signal).rolling(smoothing_window).mean()
    d_smoothed_signal = smoothed_signal.diff()
    trading = "not"
    trading_start = 0

    leverage = 0*data.copy()
    for i in range(smoothing_window,data1.shape[0]):
        leverage.iloc[i,:] = leverage.iloc[i-1,:]

        if trading=="not":

            # dynamically rebalance the signal when not trading
            fixed_beta = rolling_beta[i]
            signal = fixed_beta*data1 - data2
            smoothed_signal = pd.Series(signal).rolling(smoothing_window).mean()
            d_smoothed_signal = smoothed_signal.diff()

            if smoothed_signal[i]>0 and d_smoothed_signal[i]<0:

                leverage.iloc[i,0] = -fixed_beta / (abs(fixed_beta)+1)
                leverage.iloc[i,1] = 1 / (abs(fixed_beta)+1)

                trading = "short"
                trading_start = smoothed_signal[i]

            elif smoothed_signal[i]<0 and d_smoothed_signal[i]>0:

                fixed_beta = rolling_beta[i]
                leverage.iloc[i,0] = fixed_beta / (abs(fixed_beta)+1)
                leverage.iloc[i,1] = -1 / (abs(fixed_beta)+1)

                trading = "long"
                trading_start = smoothed_signal[i]

            else:
                leverage.iloc[i,0] = 0
                leverage.iloc[i,1] = 0


        elif trading=="long":

            # a failed trade
            if smoothed_signal[i] < trading_start:
                leverage.iloc[i,0] = 0
                leverage.iloc[i,1] = 0
                trading = "not"

            # a successful trade
            if smoothed_signal[i]>0:
                leverage.iloc[i,0] = 0
                leverage.iloc[i,1] = 0
                trading = "not"

        elif trading=="short":

            # a failed trade
            if smoothed_signal[i] > trading_start:
                leverage.iloc[i,0] = 0
                leverage.iloc[i,1] = 0
                trading = "not"

            # a successful trade
            if smoothed_signal[i]<0:
                leverage.iloc[i,0] = 0
                leverage.iloc[i,1] = 0
                trading = "not"

    return leverage

#stock1_name ="MSFT US EQUITY"
#stock2_name = "AAPL US EQUITY"
portfolioWeights = getStrategyPortfolioWeights(rolling_beta,stock1_name, stock2_name,data).fillna(0)


def backtest(pricingDF,leverageDF,start_cash):
    """Backtests pricing based on some given set of leverage. Leverage works such that it happens "overnight",
    so leverage for "today" is applied to yesterday's close price. This algo can handle NaNs in pricing data
    before a stock exists, but ffill() should be used for NaNs that occur after the stock has existed, even
    if that stock ceases to exist later."""

    pricing = pricingDF.values
    leverage = leverageDF.values

    shares = np.zeros_like(pricing)
    cash = np.zeros(pricing.shape[0])
    cash[0] = start_cash
    curr_price = np.zeros(pricing.shape[1])
    curr_price_div = np.zeros(pricing.shape[1])

    for t in range(1,pricing.shape[0]):

        if np.any(leverage[t]!=leverage[t-1]):

            # handle non-existent stock values
            curr_price[:] = pricing[t-1]     # you can multiply with this one
            curr_price[np.isnan(curr_price)] = 0
            trading_allowed = (curr_price!=0)
            curr_price_div[:] = curr_price    # you can divide with this one
            curr_price_div[~trading_allowed] = 1

            # determine new positions (warning: leverage to non-trading_allowed stocks is just lost)
            portfolio_value = (shares[t-1]*curr_price).sum()+cash[t-1]
            target_shares = trading_allowed * (portfolio_value*leverage[t]) // curr_price_div

            # rebalance
            shares[t] = target_shares
            cash[t] = cash[t-1] - ((shares[t]-shares[t-1])*curr_price).sum()

        else:

            # maintain positions
            shares[t] = shares[t-1]
            cash[t] = cash[t-1]

    returns = (shares*np.nan_to_num(pricing)).sum(axis=1)+cash
    pct_returns = (returns-start_cash)/start_cash
    return (
        pd.DataFrame( shares, index=pricingDF.index, columns=pricingDF.columns ),
        pd.Series( cash, index=pricingDF.index ),
        pd.Series( pct_returns, index=pricingDF.index)
    )

shares, cash, returns = backtest( orgin_data, portfolioWeights, 1e6 )

plt.figure(figsize = (18,8))
ax = plt.gca()
plt.title("Return Profile of Algorithm")
plt.ylabel("Percent Returns")
returns.plot(ax=ax,linewidth=3)
vals = ax.get_yticks()
ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
plt.show()
plt.savefig('2.png')


