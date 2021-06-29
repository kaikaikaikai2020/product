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
receivers='caifengsteven@gmail.com'
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
 'EURCNH', 'EURNZD', 'GBPJPY', 'GBPAUD', 'GBPCAD', 'GBPCHF', 'GBPCNH', 
 'CHFJPY', 'CHFAUD', 'CHFCAD', 'CHFNZD', 'AUDJPY', 'AUDCAD', 'AUDCNH',
 'NZDJPY', 'NZDCAD',  'CADJPY']

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
        if len(para0['sig0'])>2:
            para0['sig0'] = para0['sig0'][-10:]
            para_pool[ticker] = para0 

class memo_runnum(object):
    def __init__(self,num=0,obj=0,t0_pool=[],para_pool=[],obj_t=[]):
        self.num=num
        self.obj=0
        self.t0_pool=t0_pool
        self.para_pool=para_pool
        self.obj_t=time_use_tool()
    def update(self):
        
        self.num=self.num+1
        self.obj_t.use('%s-第%s次更新-start' % (get_now_time(),self.num))
        mark_t0 = time.time()
        #更新数据
        os.system("更新S71外汇数据.lnk")
        self.obj_t.use('%s-第%s次更新-数据更新完成' % (get_now_time(),self.num))
        #修改为多线程
        sig_info = ""
        for ticker in forex_sel:
            t0=self.t0_pool[ticker]
            para0 = self.para_pool[ticker]
            sig0_last = para0['sig0']
            x = get_forex_mintick(ticker,t0=t0,tt='2099-01-01',key_str = "tradeDate,closePrice")
            if len(x)>0:
                x.tradeDate = x.tradeDate.apply(lambda x:x[:16])
                sig0,para0= dc_bac_sig(x.tradeDate.tolist(),x.closePrice.tolist(),0.002,0,para0)
                if len(sig0)>0:
                    sig0.tradeDate0 = sig0.tradeDate0.apply(lambda x:trans_TM(x))
                    sig0['do_time'] = sig0.apply(lambda x:cal_do_time(x),axis=1)
                    para0['sig0'] = pd.concat([sig0_last,sig0])
                    para0['sig0'] = para0['sig0'][-10:]
                    #发送报告
                    tmp_str = '%0.4d-%s' % (self.num,com_sig_info(para0['sig0'],ticker))
                    if not self.obj == 0:
                        self.obj.insert(0,tmp_str)
                    else:
                        self.obj_t.use(tmp_str)
                        sig_info+=tmp_str+'\n'
                else:
                    para0['sig0'] = sig0_last
                    #tmp_str = '%0.4d-%s-keep last signal' % (self.num,ticker)
                    #if not self.obj == 0:
                    #    self.obj.insert(0,tmp_str)
                    #else:
                    #    self.obj_t.use(tmp_str)
                t0 = para0['t_end']
                t0 = slip_1M(t0)
                self.t0_pool[ticker] = t0
                self.para_pool[ticker] = para0
        if len(sig_info)>0:
            send_info(receivers,sig_info,'ZX03S71信号%s' % get_now_time())                
        print('%s-第%s次更新-complete 耗时:%s' % (get_now_time(),self.num,time.time()-mark_t0))
        return self.num
obj0 = memo_runnum(t0_pool=t0_pool,para_pool=para_pool)

def job():
    obj0.update()

scheduler = BlockingScheduler()
scheduler.add_job(job, 'interval', seconds=300)

def run_obj():
    scheduler.resume()
def stop_obj():
    scheduler.pause()

scheduler.start()
   
    

