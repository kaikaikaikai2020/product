# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 12:31:36 2021

需要做一个对接
引入多线程
多线程太复杂
@author: adair-9960
"""
'''
结果输出形式
时间 交易方向 持有时间单位
   tradeDate0 trade_dir  time_delta
0  2021-04-30         S           4

通过TK来显示

'''

#from multiprocessing.dummy import Pool as ThreadPool
from tqdm import tqdm
from yq_toolsS45 import time_use_tool
from yq_toolsS45 import read_pickle
from yq_toolsS45 import save_pickle
from bac_tool_S71 import dc_bac_sig
import pandas as pd
from yq_toolsS45 import get_forex_mintick
import datetime
import os
import time
#import tkinter
#from apscheduler.schedulers.background import BackgroundScheduler
from apscheduler.schedulers.blocking import BlockingScheduler
from M_email import send_info#(receivers,info,title_info)
receivers='nkchem09@126.com'
path0 = r'S71计算参数'
if not os.path.exists(path0):
    os.mkdir(path0)

obj_t = time_use_tool()
#时区转换
def trans_TM(x):
    x = datetime.datetime.strptime(x,'%Y-%m-%d %H:%M') \
                +datetime.timedelta(hours=8)
    return x

def slip_1M(x):
    x = datetime.datetime.strptime(x,'%Y-%m-%d %H:%M') \
                +datetime.timedelta(minutes=0.5)
    return str(x)
'''
def cal_do_time(x):
    return datetime.datetime.strptime(x.tradeDate0,'%Y-%m-%d %H:%M') \
                +datetime.timedelta(minutes=x.time_delta)
'''
def cal_do_time(x):
    return x.tradeDate0 \
                +datetime.timedelta(minutes=x.time_delta)
def com_signal(x):
    info = 'new-signal \n 计算时间: %s \n 操作方向 %s \n 持续时间:%s \n 平仓时间 %s \n\n '
    info = info % (x.tradeDate0,x.trade_dir,x.time_delta,x.do_time)
    return info
def get_now_time():
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
def com_sig_info(sig0,ticker,mod=1):
    v = sig0.iloc[-1]
    if mod==1:
        str1 = {'B':'%s-%s:平多仓一手，价格为%0.4f；买空一手，价格为：%0.4f',
                'S':'%s-%s:平空仓一手，价格为%0.4f；买多一手，价格为：%0.4f'}
        info1 = str1[v.trade_dir] % (ticker,v.tradeDate0,v.price_n,v.price_n)
    else:
        sig0['chg'] = sig0['price_n'].pct_change()
        sig0['sig'] = sig0.trade_dir.map({'B':-1,'S':1})
        sig0['chg'] = sig0['chg']*sig0.sig.shift(1)
        str1 = {'B':'%s-%s:平多仓一手，价格为%0.4f，收益为%0.4f%%：买空一手，价格为：%0.4f',
                'S':'%s-%s:平空仓一手，价格为%0.4f，收益为%0.4f%%；买多一手，价格为：%0.4f'}
        info1 = str1[v.trade_dir] % (ticker,v.tradeDate0,v.price_n,v.chg*100,v.price_n)
    print(info1)
    return info1

forex_sel = ['EURUSD', 'USDJPY', 'GBPUSD', 'AUDUSD', 'USDCAD', 'USDCHF', 'NZDUSD',
 'USDKRW', 'USDCNH', 'USDSGD', 'EURJPY', 'EURGBP', 'EURAUD', 'EURCAD', 'EURCHF',
 'EURCNH', 'EURNZD', 'GBPJPY', 'GBPAUD', 'GBPCAD', 'GBPCHF', 'GBPCNH', 'GBPNZD',
 'CHFJPY', 'CHFAUD', 'CHFCAD', 'CHFNZD', 'AUDJPY', 'AUDNZD', 'AUDCAD', 'AUDCNH',
 'NZDJPY', 'NZDCAD',  'CADJPY']

forex_sel = ['EURNZD']
#获取参数
para_pool = {}
t0_pool = {}
sig_pool = {}
for ticker in  forex_sel:
    fn = os.path.join(path0,'%s.pkl' % ticker)
    if not os.path.exists(fn):
        para0=[]
        t0 = '1990-01-01'
        sig0=pd.DataFrame()
    else:
        para0=read_pickle(fn)
        t0 = para0['t_end']
        t0 = slip_1M(t0)
        sig0=para0['sig0']
    para_pool[ticker] = para0
    t0_pool[ticker] = t0
    sig_pool[ticker] = sig0
    
#ini
obj_t.use('初始化')
os.system("更新S71外汇数据.lnk")
for ticker in tqdm(forex_sel):
    fn = os.path.join(path0,'%s.pkl' % ticker)
    t0=t0_pool[ticker]
    para0 = para_pool[ticker]
    sig0_last = sig_pool[ticker]
    x = get_forex_mintick(ticker,t0=t0,tt='2099-01-01',key_str = "tradeDate,closePrice")
    if len(x)>0:
        x.tradeDate = x.tradeDate.apply(lambda x:x[:16])
        sig0,para0= dc_bac_sig(x.tradeDate.tolist(),x.closePrice.tolist(),0.002,0,para0)
        if len(sig0)>0:
            sig0.tradeDate0 = sig0.tradeDate0.apply(lambda x:trans_TM(x))
            sig0['do_time'] = sig0.apply(lambda x:cal_do_time(x),axis=1)
            para0['sig0'] = pd.concat([sig0_last,sig0])
            para0['sig0'] = para0['sig0'][-10:]
            print(sig0.tail(10))
        else:
            para0['sig0'] = sig0_last
            print('No signal')
            
        x.tradeDate = x.tradeDate.apply(lambda x:trans_TM(x))
        print(x.tail())
        t0 = para0['t_end']
        t0 = slip_1M(t0)
        t0_pool[ticker] = t0
        para_pool[ticker] = para0   
        sig_pool[ticker] = sig0
        save_pickle(fn,para0)


    

