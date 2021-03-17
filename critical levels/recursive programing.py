# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 21:57:05 2021

@author: ASUS
"""
import pandas as pd
import numpy as np
import scipy as sp
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
import pybbg 
#import alpaca_trade_api
#service = alpaca_trade_api.REST('alpaca key', 'alpaca secret', 'https://paper-api.alpaca.markets')

bbg = pybbg.Pybbg()

fld_list = ['PX_OPEN', 'PX_HIGH', 'PX_LOW', 'PX_LAST']

def get_data(symbols, _from, _to):
    df = pd.DataFrame()
#    NY = 'America/New_York'
    start = _from
    end = _to
    #data = service.get_barset(symbols, 'day', start=start, end=end, limit=1000).df
    df = bbg.bdh(symbols, fld_list, start, end, periodselection = 'DAILY')
    #print (df)
    #df = pd.concat([df, data], axis=0)
    df['row number'] = np.arange(1, len(df)+1)
    return df
    
class TrendlineSeeker:
  
    def __init__(self, data):
        self.data = data
        self.previousTrendlineBreakingPoint = 0
        self.currentTrendlineBreakingPoint = len(data) + 1
        self.refinedTrendlines = {}
    def getCurrentRange(self):
        df = self.data.copy().loc[(self.data['row number'] >= self.previousTrendlineBreakingPoint) & (self.data['row number'] < self.currentTrendlineBreakingPoint), :]
        return df
    def getNextRange(self):
        df = self.data.loc[self.data['row number']>= self.currentTrendlineBreakingPoint, :]
        return df
    def trendlineBreakingPoints(self, currentdf, slope, intercept):
        possibleTrendBreakingPoints = currentdf.loc[(currentdf.loc[:, ("AAPL US Equity", "PX_LAST")] < 0.85*(slope*currentdf['row number'] + intercept)) | (currentdf.loc[:, ("AAPL US Equity", "PX_LAST")] > 1.15*(slope*currentdf['row number'] + intercept)), 'row number']
        return possibleTrendBreakingPoints
    def upperBoundlinesForCurrentRange(self, currentdf):
        tempdf = currentdf.copy()
        while len(tempdf) > 2:
           slope, intercept, r_value, p_value, std_err = sp.stats.linregress(x=tempdf['row number'], y=tempdf.loc[:, ("AAPL US Equity", "PX_LAST")])
           tempdf = tempdf.loc[(tempdf.loc[:, ("AAPL US Equity", "PX_LAST")]> slope * tempdf['row number'] + intercept)]
        return slope, intercept
    def refineUpperBoundlinesForCurrentRange(self, possibleUpperBoundlineBreakingPoints):
        localPossibleUpperBoundlineBreakingPoints = possibleUpperBoundlineBreakingPoints
       
        i = 1
        while len(localPossibleUpperBoundlineBreakingPoints) > 0:
            #print(len(localPossibleUpperBoundlineBreakingPoints))
            self.currentTrendlineBreakingPoint = int(localPossibleUpperBoundlineBreakingPoints[0])
            print(self.currentTrendlineBreakingPoint)
            print(self.previousTrendlineBreakingPoint)
            if self.currentTrendlineBreakingPoint - self.previousTrendlineBreakingPoint < 24:
              self.currentTrendlineBreakingPoint = len(self.data) + 1 - i
            i += 1
            currentdf = self.getCurrentRange()
            print (currentdf)
            slope, intercept = self.upperBoundlinesForCurrentRange(currentdf)
        
            localPossibleUpperBoundlineBreakingPoints = self.trendlineBreakingPoints(currentdf, slope, intercept)
            self.refinedTrendlines[str(self.previousTrendlineBreakingPoint)] = {'slope': slope, 'intercept': intercept, 'starting row': self.previousTrendlineBreakingPoint, 'ending row': self.currentTrendlineBreakingPoint - 1}
            self.previousTrendlineBreakingPoint = self.currentTrendlineBreakingPoint
            self.currentTrendlineBreakingPoint = len(self.data) + 1

    def main(self):
        i = 1
        #print('entering main')
        while True:
            currentRange = self.getCurrentRange()
            #print(currentRange)
            if len(currentRange) <= 2: break
            upperline = self.upperBoundlinesForCurrentRange(currentRange)
            #print(upperline)
            possibleUpperBoundlineBreakingPoints = self.trendlineBreakingPoints(currentRange, *upperline)
            print('possible uppper boundline breaking points')
            print(possibleUpperBoundlineBreakingPoints)
            if len(possibleUpperBoundlineBreakingPoints) == 0:
                self.refinedTrendlines[str(self.previousTrendlineBreakingPoint)] = {'slope': upperline[0], 'intercept': upperline[1], 'starting row': self.previousTrendlineBreakingPoint, 'ending row': self.currentTrendlineBreakingPoint - 1}
                break
            else: 
                print('-------')
                print(len(possibleUpperBoundlineBreakingPoints))
                print('-------')
                self.refineUpperBoundlinesForCurrentRange(possibleUpperBoundlineBreakingPoints)
                
            i += 1
            print(i)
        plt.plot(self.data.index,self.data.loc[:, ("AAPL US Equity", "PX_LAST")], label='price action')
        for key, value in self.refinedTrendlines.items():
            plt.plot(self.data.index.values[value['starting row']:value['ending row']], value['slope']*self.data['row number'][value['starting row']:value['ending row']] + value['intercept'])
            plt.legend(loc='best')
            plt.show()
if __name__ == '__main__':
   #data = get_data('AAPL', '2019-08-29', '2021-02-21')
   data = get_data(['AAPL US Equity'], '20190829', '20210221')
   #print(data)
   #print(len(data))
   ts = TrendlineSeeker(data)
   ts.main()