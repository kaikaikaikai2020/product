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
   
db1 = create_db('polygon_stock_minute')

fn1 = 'SPY'
fn2='QQQ'

sql_str= 'select tradeDate as date,closePrice as close from %s order by tradeDate'

x=pd.read_sql(sql_str % fn1,db1)
tmp = x.date.apply(lambda x:True if np.mod(x.minute,30)==0 else False)
x=x[tmp]
x.date = x.date.astype(str)

y=pd.read_sql(sql_str % fn2,db1)
tmp = y.date.apply(lambda x:True if np.mod(x.minute,30)==0 else False)
y=y[tmp]
y.date = y.date.astype(str)
do_cal([fn1,fn2,x,y,None,False])




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