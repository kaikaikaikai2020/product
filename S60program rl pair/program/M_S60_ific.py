# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 11:17:58 2020

@author: adair-9960
"""

import pandas as pd
import numpy as np
from S60_bac_tool import do_cal
#from S60_bac_tool import multiprocessing,num_core
from yq_toolsS45 import create_db
train_len=None
   
db1 = create_db('future_min_data')
db2 = create_db('polygon')

#fn1_pool = ['SPY','ES','ES','ES_day']
#fn2_pool = ['QQQ','VX','NQ','VX_day']
fn1_pool = ['ccfx_ic']
fn2_pool = ['ccfx_if']

sql_str1= 'select tradingdate as date,close from %s order by tradingdate'
sql_str2 = 'select tradeDate as date,closePrice as close from usastock_day where ticker = "%s" order by tradeDate'

if __name__=='__main__':
    f = []
    for fn1,fn2 in zip(fn1_pool,fn2_pool):
        if 'day' in fn1:
            sub_fn1 = fn1.split('_day')[0]
            sub_fn2 = fn1.split('_day')[0]
            sql_str = sql_str2
            x=pd.read_sql(sql_str % sub_fn1,db2)
            x.date = x.date.astype(str)
            y=pd.read_sql(sql_str % sub_fn2,db2)
            y.date = y.date.astype(str)
        else:
            sql_str = sql_str1
            x=pd.read_sql(sql_str % fn1,db1)
            tmp = x.date.apply(lambda x:True if np.mod(x.minute,30)==0 else False)
            x=x[tmp]
            x.date = x.date.astype(str)
            
            y=pd.read_sql(sql_str % fn2,db1)
            tmp = y.date.apply(lambda x:True if np.mod(x.minute,30)==0 else False)
            y=y[tmp]
            y.date = y.date.astype(str)
        f.append([fn1,fn2,x,y,None,True])
    do_cal(f[0])





"""




if __name__=='__main__':
    f = []
    data = updateS60data()
    for info in data.keys():
        fn1,fn2=info.split('-')
        x,y=data[info]
        f.append([fn1,fn2,x,y,train_len,False])
    pool = multiprocessing.Pool(min(num_core,len(f))) 
    pool.map(do_cal, f)
    pool.close()
    pool.join()
"""