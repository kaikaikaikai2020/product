# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 11:17:58 2020

@author: adair-9960
"""

import pandas as pd
import numpy as np
from S60_bac_tool import do_cal
from S60_bac_tool import multiprocessing,num_core
from yq_toolsS45 import create_db
train_len=None
   
db1 = create_db('us_stock')

fn1_pool = ['ES','ES']
fn2_pool=['NQ','VX']

sql_str= """select concat(UTCDate, ' ', UTCTimeBarStart) as date,closePrice as close 
from futureminute_m where Ticker  like '%s%%' order by date;"""

if __name__ == '__main__':
    f = []
    for fn1,fn2 in zip(fn1_pool,fn2_pool):
        x=pd.read_sql(sql_str % fn1,db1)
        x.date = pd.to_datetime(x.date)
        tmp = x.date.apply(lambda x:True if np.mod(x.minute,30)==0 else False)
        x=x[tmp]
        x.date = x.date.astype(str)
        
        y=pd.read_sql(sql_str % fn2,db1)
        y.date = pd.to_datetime(y.date)
        tmp = y.date.apply(lambda x:True if np.mod(x.minute,30)==0 else False)
        y=y[tmp]
        y.date = y.date.astype(str)
        f.append([fn1,fn2,x,y,None,True])
    pool = multiprocessing.Pool(min(num_core,len(f))) 
    pool.map(do_cal, f)
    pool.close()
    pool.join()




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