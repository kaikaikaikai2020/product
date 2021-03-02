# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 11:17:58 2020

@author: adair-9960
"""

from S60_bac_tool import updateS60data,do_cal
from S60_bac_tool import multiprocessing,num_core
train_len=800
   

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