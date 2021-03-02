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

fn1_pool = ['ES_day','ES_day']
fn2_pool=['NQ_day','VX_day']

if __name__ == '__main__':
    f = []
    for fn1,fn2 in zip(fn1_pool,fn2_pool):
        
        sub_fn1 = fn1.split('_day')[0]
        
        x=pd.read_csv('%s.csv' % sub_fn1,dtype={'date':str})
        sub_fn2 = fn2.split('_day')[0]
        y=pd.read_csv('%s.csv' % sub_fn2,dtype={'date':str})
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