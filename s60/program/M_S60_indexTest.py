# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 11:17:58 2020

@author: adair-9960
"""

import pandas as pd
from yq_toolsS45 import engine,list_intersec
from S60_bac_tool_TMP import do_cal
from S60_bac_tool_TMP import multiprocessing,num_core
train_len=800
   

if __name__=='__main__':
    data0=pd.read_csv('T000300.csv',dtype={'ticker':int,'ticker2':int,'t':str})
    data0.ticker1 = data0.ticker1.apply(lambda x:'%0.6d' % x)
    data0.ticker2 = data0.ticker2.apply(lambda x:'%0.6d' % x)
    fn1 = data0.ticker1.tolist()
    fn2 = data0.ticker2.tolist()
    tpool=data0.t.tolist()
    f = []
    sql_str = 'select tradeDate as date,closePrice as close from yq_MktEqudAdjAfGet where ticker="%s" order by tradeDate'
    for sub_fn1,sub_fn2,tt in zip(fn1,fn2,tpool):
        x=pd.read_sql(sql_str % sub_fn1,engine)
        x.date = x.date.astype(str)
        
        y=pd.read_sql(sql_str % sub_fn2,engine)
        y.date = y.date.astype(str)
    
        tmp = list_intersec(x.date.tolist(),y.date.tolist())
        train_len = tmp.index(tt)
        
        if (train_len)>400:
            f.append([sub_fn1,sub_fn2,x,y,tt,True])
            
    pool = multiprocessing.Pool(min(num_core,len(f))) 
    pool.map(do_cal, f)
    pool.close()
    pool.join()