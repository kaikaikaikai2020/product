# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 08:37:24 2020

@author: adair-9960
"""

import pandas as pd
from yq_toolsS45 import time_use_tool

def updateS60data():
    obj_clock=time_use_tool()
    def trans_sub_x(sub_x):
        sub_x = sub_x.copy()
        sub_x.dropna(how='all',inplace=True)
        sub_x.dropna(how='all',axis='columns',inplace=True)
        sub_x=sub_x.iloc[2:]
        sub_x.columns=['date','close']
        sub_x['date'] = sub_x.date.apply(lambda x:x.date()).astype(str)
        return sub_x
    
    x = pd.read_excel('pair trading data CN HK KR TW.xlsx')
    x.dropna(how='all',inplace=True)
    x.dropna(how='all',axis='columns',inplace=True)
    x=x.iloc[2:]
    x.iloc[0].fillna(method='bfill',inplace=True)
    x.reset_index(drop =True,inplace=True)
    x.drop(1,inplace=True)
    x=x.T
     
    data = {}   
    while len(x)>0:
        sub_x = x.iloc[0:4].T
        x=x.iloc[4:]
        
        sub_x1=sub_x[sub_x.columns[0:2]]
        sub_x1 = trans_sub_x(sub_x1)
        sub_x2=sub_x[sub_x.columns[2:4]]
        sub_x2 = trans_sub_x(sub_x2)
        sub_ticker=sub_x.iloc[0].unique().tolist()
        sub_ticker = [i.split(' ')[0] for i in sub_ticker]
        data['-'.join(sub_ticker)] = [sub_x1,sub_x2]
    obj_clock.use('Time used')
    return data