# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 14:35:35 2020

@author: adair-9960
"""
import pandas as pd
import os
pn0 = r'HA_result'
from yq_toolsS45 import read_pickle

if not os.path.exists(pn0):
    os.mkdir(pn0)
if __name__ == '__main__':
    x=pd.read_csv('stockListHA.csv',dtype={'f1':str,'f2':str})
    fn_a = zip(x.f1.tolist(),x.f2.tolist())
    
    for sub_fn in fn_a:
        fn1,fn2 = sub_fn
        com_fn1 =os.path.join(pn0, 'bac-%s-%s.pkl' % (fn1,fn2))
        pnl = read_pickle(com_fn1)
        com_fn1 =os.path.join(pn0, 'bac-%s-%s.csv' % (fn1,fn2))
        pnl.to_csv(com_fn1)