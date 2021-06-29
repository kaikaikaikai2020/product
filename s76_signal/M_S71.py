# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 15:24:52 2021
国内指数 
us_forex 十万分之五
美股 千分之1.2
@author: adair-9960
"""
import matplotlib.pyplot as plt
from bac_tool_S71 import dc_bac_v1
from yq_toolsS45 import engine
import pandas as pd
from yq_toolsS45 import create_db
from tqdm import tqdm
from yq_toolsS45 import save_pickle,read_pickle
import os

save_dir = 'S71_cal_result'
if not os.path.exists(save_dir):
    os.mkdir(save_dir)


eg_tdx=create_db('pytdx_data')
eg_28 = create_db('s28')
eg_csi_f = create_db('future_min_data') #csi_future_min
eg_us_minute = create_db('polygon_stock_minute')
eg_us_forex = create_db('polygon_fx_minute')
eg_is_index = create_db('s40_america')
eg_us_future = create_db('us_stock')
eg_polygon = create_db('polygon')

#tn yq_MktMFutdGet  期货主力合约
#tn 


fee={'usstock':1.2/1000,'csi_index':5/10000,'csi_future':3/10000,
     'us_forex':5/100000,'is_index':5/10000,'us_future':5/100000,
     'us_forex_day':5/100000}


def rep_day_inf(t2):
    t2 = t2.pct_change()
    t2['t1'] = pd.to_datetime(t2.index.tolist())
    t2['t1'] = t2.t1.apply(lambda x:x.date())
    t2['t2'] = t2.t1.shift(1)
    c = ['dc', 'dc_close', 'bac_dc_dc_close']
    t2[t2.t1!=t2.t2] = 0
    t2 = t2[c]
    return t2

#外汇
def us_forex():
    fn = os.path.join(save_dir,'us_forex.pkl')
    if os.path.exists(fn):
       Y = read_pickle(fn)
    else:
        Y = []
        for code in tqdm(['eurusd','eurjpy','eurgbp','euraud','eurchf','eurcad','gbpusd',
                 'gbpchf','gbpjpy','gbpcad','gbpaud','usdjpy','audusd','nzdusd', 
                 'audcad','audjpy','audnzd','nzdjpy','chfjpy','usdchf','usdcad',
                 'usdcnh','usdmxn','usdtry','usdkrw']):
            sql_str_select_data1 = '''select tradeDate,closePrice from %s order by tradeDate'''
            sql_str_select_data1 = sql_str_select_data1 % (code)
            x = pd.read_sql(sql_str_select_data1, eg_us_forex)
            if len(x)>100:
                tmp = dc_bac_v1(x.tradeDate.tolist(),x.closePrice.tolist(),0.002,0)[0]
                Y.append({code:tmp})
                save_pickle(fn,Y)
    X = {}
    for i in Y:
        X.update(i)
    for t1 in X.keys():
        t2 = X[t1].copy()
        t2 = rep_day_inf(t2)
        t2 = ((t2-fee['us_forex']).fillna(0)+1).cumprod()
        t2.plot(figsize=(10,5),title=t1)
        plt.show()
        
#国内指数
def P1_csi_index():
    for sub_index in [300,500,50]:
        sql_tmp = 'select tradingdate,closeprice as close from wind_%d_update order by tradingdate' % sub_index
        x = [pd.read_sql(sql_tmp,eg_28)]
        
        
        sql_tmp = 'select tradingdate,close from tdx_min_%d where tradingdate >"%s" order by tradingdate'
        x.append(pd.read_sql(sql_tmp % (sub_index,str(x[0].tradingdate.max())),eg_tdx))
        x = pd.concat(x)
        y = dc_bac_v1(x.tradingdate.tolist(),x.close.tolist(),0.05)[0]
        y = ((y.pct_change()-fee['csi_index']).fillna(0)+1).cumprod()
        y.plot(figsize=(10,5),title=("csi%d" % sub_index))
        #return y

#国内期货
def P1_csi_future():
    fn = os.path.join(save_dir,'csi_future_0.004_1.pkl')
    if os.path.exists(fn):
        Y = read_pickle(fn)
    else:
        tn=pd.read_sql('show tables',eg_csi_f)
        tn=tn[tn.columns[0]].tolist()
        tn=[i for i in tn if i[-2:]!='dt']
        sql_f1 = 'select tradingdate,close from %s order by tradingdate'
        Y = []
        for sub_tn in tqdm(tn):
            x = pd.read_sql(sql_f1 % sub_tn,eg_csi_f)
            if len(x)>100:
                tmp = dc_bac_v1(x.tradingdate.tolist(),x.close.tolist(),0.004,0)[0]
                Y.append({sub_tn:tmp})
        save_pickle(fn,Y)
    X = {}
    for i in Y:
        X.update(i)
    for t1 in X.keys():
        t2 = X[t1].copy()
        t2 = rep_day_inf(t2)
        #t2 = ((t2-fee['csi_future']).fillna(0)+1).cumprod()        
        #t2 = t2.pct_change().fillna(0).cumsum()
        t2 = (t2-fee['csi_future']).fillna(0).cumsum()
        if t2.dc.max()>1:
            t2.plot(figsize=(10,5),title=t1)
            plt.show()
        
#美国股票分钟
def us_stock():
    for ticker in ['aapl','qqq']:
        sql_str = 'select tradeDate,closePrice from %s order by tradeDate' % ticker
        x=pd.read_sql(sql_str,eg_us_minute)
        y = dc_bac_v1(x.tradeDate.tolist(),x.closePrice.tolist(),0.03)[0]    
        y = rep_day_inf(y)
        y = ((y-fee['usstock']).fillna(0)+1).cumprod()
        y.plot(figsize=(10,5),title=ticker)    

#国际指数
def is_index():
    fn = os.path.join(save_dir,'is_index千2.pkl')
    if os.path.exists(fn):
       Y = read_pickle(fn)
    else:
        Y = []
        tns = pd.read_sql('show tables',eg_is_index)
        tns = tns[tns.columns[0]].tolist()
        for code in tqdm(tns):
            sql_str_select_data1 = '''select tradingdate as tradeDate,closePrice from %s order by tradingdate'''
            sql_str_select_data1 = sql_str_select_data1 % (code)
            x = pd.read_sql(sql_str_select_data1, eg_is_index)
            if len(x)>100:
                tmp = dc_bac_v1(x.tradeDate.tolist(),x.closePrice.tolist(),0.01,0)[0]
                Y.append({code:tmp})
        save_pickle(fn,Y)
    X = {}
    for i in Y:
        X.update(i)
    for t1 in X.keys():
        t2 = X[t1].copy()
        t2 = t2.pct_change()
        t2 = ((t2-fee['is_index']).fillna(0)+1).cumprod()
        #t2 = (t2.fillna(0)+1).cumprod()
        t2.plot(figsize=(10,5),title=t1)
        plt.show()  
#美国期货
def us_future4():
    fn = os.path.join(save_dir,'us_future.pkl')
    if os.path.exists(fn):
       Y = read_pickle(fn)
    else:        
        tickers = ['ES','NQ','ZB','VX']
        sql_tmp = 'select ticker,CONCAT(UTCDate," ",UTCTimeBarStart) as tradeDate,closePrice from futureminute_m where ticker like "%s%%"  order by tradeDate'
        Y = []
        for ticker in tickers:
            x = pd.read_sql(sql_tmp % ticker,eg_us_future)
            if len(x)>100:
                tmp = dc_bac_v1(x.tradeDate.tolist(),x.closePrice.tolist(),0.05,0)[0]
                tmp = tmp.pct_change()
                tmp = tmp.merge(x.set_index('tradeDate'),left_index=True,right_index=True)
                var = ['dc', 'dc_close', 'bac_dc_dc_close']
                a = tmp[var].copy()
                a[tmp.ticker!=tmp.ticker.shift(1)] = 0
                tmp = a.copy()
                Y.append({ticker:tmp})
        save_pickle(fn,Y)
    X = {}
    for i in Y:
        X.update(i)
    for t1 in X.keys():
        t2 = X[t1].copy()
        #t2 = rep_day_inf(t2)
        #t2 = t2.pct_change()
        t2 = ((t2-fee['us_future']).fillna(0)+1).cumprod()
        #t2 = (t2.fillna(0)+1).cumprod()
        t2.index = t2.index.astype(str)
        t2.index=pd.to_datetime(t2.index.tolist())
        t2.plot(figsize=(10,5),title=t1)
        plt.show() 
        
#国内日度频率数据
def us_forex_day():
    fn = os.path.join(save_dir,'us_forex_d.pkl')
    if os.path.exists(fn):
       Y = read_pickle(fn)
    else:
        Y = []
        for code in tqdm(['eurusd','eurjpy','eurgbp','euraud','eurchf','eurcad','gbpusd',
                 'gbpchf','gbpjpy','gbpcad','gbpaud','usdjpy','audusd','nzdusd', 
                 'audcad','audjpy','audnzd','nzdjpy','chfjpy','usdchf','usdcad',
                 'usdcnh','usdmxn','usdtry','usdkrw']):
            sql_str_select_data1 = '''select tradeDate,closePrice from forex_day where ticker =  "%s" order by tradeDate'''
            sql_str_select_data1 = sql_str_select_data1 % (code)
            x = pd.read_sql(sql_str_select_data1, eg_polygon)
            if len(x)>100:
                tmp = dc_bac_v1(x.tradeDate.tolist(),x.closePrice.tolist(),0.002,0)[0]
                Y.append({code:tmp})
                save_pickle(fn,Y)
    X = {}
    for i in Y:
        X.update(i)
    for t1 in X.keys():
        t2 = X[t1].copy()
        #t2 = rep_day_inf(t2)
        t2 = ((t2.pct_change()-fee['us_forex_day']).fillna(0)+1).cumprod()
        t2.plot(figsize=(10,5),title=t1)
        plt.show()
        
def P1_csi_index_d():
    for sub_index in ['000905','000300','000016']:
        sql_tmp = 'select tradeDate,closeIndex as close from yq_index where symbol ="%s" order by tradeDate' % sub_index
        x = pd.read_sql(sql_tmp,engine)
        y = dc_bac_v1(x.tradeDate.tolist(),x.close.tolist(),0.05)[0]
        y = ((y.pct_change()-fee['csi_index']).fillna(0)+1).cumprod()
        #y = ((y.pct_change()-0).fillna(0)+1).cumprod()
        y.plot(figsize=(10,5),title=("csi%s" % sub_index))
        
def P1_csi_future_day():    
    fn = os.path.join(save_dir,'csi_future_d2.pkl')
    if os.path.exists(fn):
        Y = read_pickle(fn)
    else:
        ticker ='select distinct(concat(exchangeCD,"-",contractObject)) from yq_FutuGet'
        ticker = pd.read_sql(ticker,engine)
        ticker = ticker[ticker.columns[0]].tolist()
        
        sql_f1 = '''select ticker,tradeDate as tradingdate,closePrice as close from yq_MktMFutdGet 
                where exchangeCD = "%s" and contractObject="%s" and mainCon = 1 order by tradingdate'''
        Y = []
        #ticker = ['CCFX-IC','CCFX-IH','CCFX-IF']
        for sub_tn in tqdm(ticker):
            v1,v2=sub_tn.split('-')
            x = pd.read_sql(sql_f1 % (v1,v2),engine)
            if len(x)>100:
                tmp = dc_bac_v1(x.tradingdate.tolist(),x.close.tolist(),0.004,1)[0]
                #换约修正
                tmp = tmp.pct_change()
                tmp = tmp.merge(x.set_index('tradingdate'),left_index=True,right_index=True)
                var = ['dc', 'dc_close', 'bac_dc_dc_close']
                a = tmp[var].copy()
                a[tmp.ticker!=tmp.ticker.shift(1)] = 0
                tmp = a.copy()
                Y.append({sub_tn:tmp})
        save_pickle(fn,Y)
    X = {}
    for i in Y:
        X.update(i)
    for t1 in X.keys():
        t2 = X[t1].copy()
        t2 = (t2-fee['csi_future']).fillna(0).cumsum()
        t2.plot(figsize=(10,5),title=t1)
        plt.show()
        
if __name__ == "__main__":
    #外汇
    us_forex()
    #国内指数
    #国内指数
    P1_csi_index()
    #国内期货
    P1_csi_future()
    #美国股票分钟
    us_stock()
    #国际指数
    is_index()
    #美国期货
    us_future4()
    #日度频率数据
    #外汇
    us_forex_day()
    #指数
    P1_csi_index_d()
    #期货
    P1_csi_future_day()