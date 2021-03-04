# -*- coding: utf-8 -*-
"""
Created on Tue May 26 16:45:10 2020

@author: adair2019
"""

import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sqlalchemy import create_engine
import json
from datetime import date,datetime
import pymysql
import warnings
import sys
import time
import os
import multiprocessing
from multiprocessing.dummy import Pool as ThreadPool

from time import strftime, localtime
from datetime import timedelta
import calendar

import pickle
warnings.filterwarnings('ignore')
num_core = multiprocessing.cpu_count()
#must be set before using
with open('para.json','r',encoding='utf-8') as f:
    para = json.load(f)
    
pn = para['yuqerdata_dir']

user_name = para['mysql_para']['user_name']
pass_wd = para['mysql_para']['pass_wd']
port = para['mysql_para']['port']
host='localhost'

db_name1 = 'yuqerdata'
#eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name)
eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name1)
engine = create_engine(eng_str,pool_recycle=10600, pool_size=100, max_overflow=20)

db_tdx = 'pytdx_data'
#eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name)
eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_tdx)
eg_tdx = create_engine(eng_str)

db_name23 = 'S23'
#eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name)
eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name23)
eg_23 = create_engine(eng_str)

db_name31 = 'S31'
#eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name)
eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name31)
eg_31 = create_engine(eng_str)

db_name40 = 'S40'
#eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name)
eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name40)
eg_40 = create_engine(eng_str)
#
db_40index = 's40_america'
#eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name)
eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_40index)
eg_40index = create_engine(eng_str)

db_name42 = 'S42'
#eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name)
eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name42)
engine42 = create_engine(eng_str)

db_name46 = 'S46'
#eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name)
eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name46)
engine46 = create_engine(eng_str)

db_name48 = 'S48'
#eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name)
eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name48)
engine48 = create_engine(eng_str)

db_name49 = 'S49'
eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name49)
engine49 = create_engine(eng_str)

db_name37 = 'S37'
eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name37)
engine37 = create_engine(eng_str)

db_name_us = 'us_stock'
#eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name)
eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name_us)
engine_us = create_engine(eng_str)

db_name_polygon_fx_minute= 'polygon_fx_minute'
eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name_polygon_fx_minute)
engine_db_name_polygon_fx_minute = create_engine(eng_str)

#polygon
db_name_polygon= 'polygon'
eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name_polygon)
engine_db_polygon = create_engine(eng_str)

db_name_yq_cub = 'yuqer_cubdata_update'
eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name_yq_cub)
engine_yq_cub = create_engine(eng_str)

db_para = 'parapool'
eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_para)
engine_para = create_engine(eng_str)

db_usindexmin = 'foreign_index_min'
eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_usindexmin)
eb_usindminu = create_engine(eng_str)

db_ak = 'aksharedata'
#eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name)
eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_ak)
eg_ak = create_engine(eng_str)

def create_db(db_ak):
    #eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_name)
    eng_str='mysql+pymysql://%s:%s@localhost:%d/%s?charset=utf8' % (user_name,pass_wd,port,db_ak)
    eg_ak = create_engine(eng_str)
    return eg_ak


sql_str_select_data1 = '''select %s from yq_dayprice where symbol="%s" and tradeDate>="%s"
    and tradeDate<="%s" order by tradeDate'''
sql_str_select_data2 = '''select %s from MktEqudAdjAfGet where ticker="%s" and tradeDate>="%s"
    and tradeDate<="%s" order by tradeDate'''

def get_inidata(tn,key_str='tradeDate',eg=engine):
    sql_str = 'select %s from %s order by %s desc limit 1'
    t = pd.read_sql(sql_str % (key_str,tn,key_str),eg)
    if len(t)>0:
        t = t[t.columns[0]].astype(str).values[0]
    else:
         t = '1989-01-01'
    return t

def get_file_name(file_dir,file_type):
    L=[]
    L_s = []   
    for root, dirs, files in os.walk(file_dir):  
        for file in files:  
            if os.path.splitext(file)[1] == file_type:  
                L.append(os.path.join(root, file))  
                L_s.append(file)
    return L,L_s
#列表转带引号str
def list_to_str_f(week_end_list):
    return  '''"'''+'","'.join(week_end_list)+'''"'''
#创建表格并分区
def create_table_update(db_name,tn_name,var_name,var_type,key_str,p_num=1,partions_str='ticker'):
    #连接本地数据库
    db = pymysql.connect(host,user_name,pass_wd,db_name)
    #创建游标
    cursor = db.cursor()
    #创建
    var_info=''
    for id,sub_var in enumerate(var_name):
        var_info=var_info + sub_var + ' ' + var_type[id] + ','
    var_info = var_info[:-1]    
    if len(key_str)>0:
        sql = 'create table  `%s`(%s,primary key(%s)) partition by key(%s) partitions %d' % (tn_name,var_info,key_str,partions_str,p_num)    
    else:
        sql = 'create table  `%s`(%s)' % (tn_name,var_info)  

    try:
        # 执行SQL语句
        cursor.execute(sql)
        print("创建数据库成功")
    except Exception as e:
        print("创建数据库失败：case%s"%e)
    finally:
        #关闭游标连接
        cursor.close()
        # 关闭数据库连接
        db.close()
        
        
def create_table(db_name,tn_name,var_name,var_type,key_str):
    #连接本地数据库
    db = pymysql.connect(host,user_name,pass_wd,db_name)
    #创建游标
    cursor = db.cursor()
    #创建
    var_info=''
    for id,sub_var in enumerate(var_name):
        var_info=var_info + sub_var + ' ' + var_type[id] + ','
    var_info = var_info[:-1]    
    if len(key_str)>0:
        sql = 'create table  `%s`(%s,primary key(%s))' % (tn_name,var_info,key_str)    
    else:
        sql = 'create table  `%s`(%s)' % (tn_name,var_info)  

    try:
        # 执行SQL语句
        cursor.execute(sql)
        print("创建数据库成功")
    except Exception as e:
        print("创建数据库失败：case%s"%e)
    finally:
        #关闭游标连接
        cursor.close()
        # 关闭数据库连接
        db.close()

def do_sql_order(order_str,db_name):
    db = pymysql.connect(host,user_name,pass_wd,db_name)
    #创建游标
    cursor = db.cursor()
    try:
        # 执行SQL语句
        cursor.execute(order_str)
        print("执行mysql命令成功")
    except Exception as e:
        print("执行mysql命令失败：case%s"%e)
    finally:
        #关闭游标连接
        cursor.close()
        # 关闭数据库连接
        db.close()
        
        
def get_a_stock_tradeDate_S45(index,begin='2000-01-01',end='2033-01-01'):
    #每日数据
    sql_str_index = '''select symbol,tradeDate,closePrice  as closeIndex,turnoverVol
    from yq_dayprice     where symbol = "%s"  and tradeDate>="%s" and 
    tradeDate<="%s" order by tradeDate'''
    sql_str_index = sql_str_index % (index,begin,end)
    hs300_index = pd.read_sql(sql_str_index,engine)
    hs300_index = hs300_index.sort_values('tradeDate')
    
    #后复权系数 MktEqudAdjAfGet
    sql_str_fq = """select tradeDate,accumAdjFactor from MktEqudAdjAfGet where 
    ticker = "%s" order by tradeDate"""
    sql_str_fq = sql_str_fq % index
    y = pd.read_sql(sql_str_fq,engine)
    y=pd.merge(hs300_index,y,on=['tradeDate'])
    y['closeIndex'] = y['closeIndex']*y['accumAdjFactor']
    y.drop(columns=['accumAdjFactor'],inplace=True)
    return y        
# 指数数据数据的起始与终止时间
def get_index_tradeDate(index,begin,end,var_name = '*'):
    if isinstance(var_name,str):
        var_sel = var_name
    else:
        var_sel = ','.join(var_name)
            
    sql_str_index = '''select %s from yq_index where symbol = "%s" and tradeDate>="%s" 
        and tradeDate<="%s" order by tradeDate'''
    sql_str_index = sql_str_index % (var_sel,index,begin,end)
    hs300_index = pd.read_sql(sql_str_index,engine)
    hs300_index = hs300_index.sort_values('tradeDate')
    return hs300_index
#指数日度数据
def MktIdxdGet(index,begin,end,var_name = '*'):
    if isinstance(var_name,str):
        var_sel = var_name
    else:
        var_sel = ','.join(var_name)
            
    sql_str_index = '''select %s from yq_index where symbol = "%s" and tradeDate>="%s" 
        and tradeDate<="%s" order by tradeDate'''
    sql_str_index = sql_str_index % (var_sel,index,begin,end)
    hs300_index = pd.read_sql(sql_str_index,engine)
    hs300_index = hs300_index.sort_values('tradeDate')
    return hs300_index
#指数日度数据 多个
def MktIdxdGetCom(indexlist,begin,end,var_name = '*'):
    if isinstance(indexlist,list):
        index=list_to_str_f(indexlist)
    else:
        index = indexlist
    if isinstance(var_name,str):
        var_sel = var_name
    else:
        var_sel = ','.join(var_name)
            
    sql_str_index = '''select %s from yq_index where symbol in (%s) and tradeDate>="%s" 
        and tradeDate<="%s" order by tradeDate'''
    sql_str_index = sql_str_index % (var_sel,index,begin,end)
    hs300_index = pd.read_sql(sql_str_index,engine)
    hs300_index = hs300_index.sort_values('tradeDate')
    return hs300_index
#指数月度数据
def MktIdxmGet(index,begin,end,var_name = '*'):
    if isinstance(var_name,str):
        var_sel = var_name
    else:
        var_sel = ','.join(var_name)
            
    sql_str_index = '''select %s from yq_index_month where symbol = "%s" and endDate>="%s" 
        and endDate<="%s" order by endDate'''
    sql_str_index = sql_str_index % (var_sel,index,begin,end)
    hs300_index = pd.read_sql(sql_str_index,engine)
    hs300_index = hs300_index.sort_values('endDate')
    return hs300_index
    
def get_cf_future_update(index,begin='2000-01-01',end='2033-01-01',var_name = '*'):
    if isinstance(var_name,str):
        var_sel = var_name
    else:
        var_sel = ','.join(var_name)
    sql_str_index = '''select %s from yq_MktMFutdGet 
    where contractObject = "%s" and mainCon = 1 and tradeDate>="%s" and 
    tradeDate<="%s" order by tradeDate'''
    sql_str_index = sql_str_index % (var_sel,index,begin,end)
    hs300_index = pd.read_sql(sql_str_index,engine)
    hs300_index = hs300_index.sort_values('tradeDate')
    return hs300_index

def get_cf_future_tradeDate(index,begin='2000-01-01',end='2033-01-01'):
    sql_str_index = '''select contractObject as symbol,tradeDate,openPrice as openIndex,highestPrice as highestIndex,
    lowestPrice as lowestIndex,closePrice  as closeIndex,turnoverVol,chgPct from yq_MktMFutdGet 
    where contractObject = "%s" and mainCon = 1 and tradeDate>="%s" and 
    tradeDate<="%s" order by tradeDate'''
    sql_str_index = sql_str_index % (index,begin,end)
    hs300_index = pd.read_sql(sql_str_index,engine)
    hs300_index = hs300_index.sort_values('tradeDate')
    return hs300_index

def get_a_stock_tradeDate(index,begin='2000-01-01',end='2033-01-01'):
    #每日数据
    sql_str_index = '''select symbol,tradeDate,openPrice as openIndex,highestPrice as highestIndex,
    lowestPrice as lowestIndex,closePrice  as closeIndex from yq_dayprice
    where symbol = "%s"  and tradeDate>="%s" and 
    tradeDate<="%s" order by tradeDate'''
    sql_str_index = sql_str_index % (index,begin,end)
    hs300_index = pd.read_sql(sql_str_index,engine)
    hs300_index = hs300_index.sort_values('tradeDate')
    
    #后复权系数 MktEqudAdjAfGet
    sql_str_fq = """select tradeDate,accumAdjFactor from MktEqudAdjAfGet where 
    ticker = "%s" order by tradeDate"""
    sql_str_fq = sql_str_fq % index
    y = pd.read_sql(sql_str_fq,engine)
    y=pd.merge(hs300_index,y,on=['tradeDate'])
    y['openIndex'] = y['openIndex']*y['accumAdjFactor']
    y['highestIndex'] = y['highestIndex']*y['accumAdjFactor']
    y['lowestIndex'] = y['lowestIndex']*y['accumAdjFactor']
    y['closeIndex'] = y['closeIndex']*y['accumAdjFactor']
    return y

def get_exchange_tradeDate(index,begin='2000-01-01',end='2033-01-01'):
    #每日数据
    sql_str_index = '''select symbol,tradingdate as tradeDate,openPrice as openIndex,highestPrice as highestIndex,
    lowestPrice as lowestIndex,closePrice  as closeIndex,turnoverVol from exchange_dayly
    where symbol = "%s"  and tradingdate>="%s" and 
    tradingdate<="%s" order by tradingdate'''
    sql_str_index = sql_str_index % (index,begin,end)
    hs300_index = pd.read_sql(sql_str_index,engine42)
    hs300_index = hs300_index.sort_values('tradeDate')    
    return hs300_index

def get_exchange_tradeDate_update(index,begin='2000-01-01',end='2033-01-01',var_name = '*'):
    #每日数据
    sql_str_index = '''select %s from exchange_dayly
    where symbol = "%s"  and tradingdate>="%s" and 
    tradingdate<="%s" order by tradingdate'''
    sql_str_index = sql_str_index % (var_name,index,begin,end)
    hs300_index = pd.read_sql(sql_str_index,engine42)
    hs300_index = hs300_index.sort_values('tradeDate')    
    return hs300_index

#dowjones data
def get_dowjones_tradeDate(index,begin='2000-01-01',end='2033-01-01'):
    #每日数据
    sql_str_index = '''select symbol,tradeDate,openPrice as openIndex,highestPrice as highestIndex,
    lowestPrice as lowestIndex,closePrice  as closeIndex, totalVolume as turnoverVol from dowjones_dayly
    where symbol = "%s"  and tradeDate>="%s" and 
    tradeDate<="%s" order by tradeDate'''
    sql_str_index = sql_str_index % (index,begin,end)
    hs300_index = pd.read_sql(sql_str_index,engine42)
    hs300_index = hs300_index.sort_values('tradeDate')    
    return hs300_index

#美股后复权数据
def get_american_stock_tradeDate(index,begin='2000-01-01',end='2033-01-01'):
    #每日数据
    sql_str_index = '''select symbol,tradingdate as tradeDate,openprice_adj as openIndex,highprice_adj as highestIndex,
    lowprice_adj as lowestIndex,closeprice_adj  as closeIndex,volume_adj as turnoverVol from us_stock_daytick
    where symbol = "%s"  and tradingdate>="%s" and 
    tradingdate<="%s" order by tradingdate'''
    sql_str_index = sql_str_index % (index,begin,end)
    hs300_index = pd.read_sql(sql_str_index,engine_us)
    return hs300_index
#获取成分股
def get_IdxCons(intoDate,ticker='000300'):
    #nearst 时间
    sql_str1 = '''select symbol from yuqerdata.IdxCloseWeightGet where ticker = "%s"
            and tradingdate = (select tradingdate from yuqerdata.IdxCloseWeightGet where 
        ticker="%s" and tradingdate<="%s"  order by tradingdate desc limit 1)''' %(ticker,
        ticker,intoDate)
    x = pd.read_sql(sql_str1,engine)
    x = x['symbol'].values   
    return x

#日线数据
def chg_factor(ticker = '000005',begin = '2001-01-01' ,end = '2090-01-01' , 
               field = [u'symbol',  u'tradeDate', u'openPrice',
                        u'highestPrice', u'lowestPrice', u'closePrice', u'turnoverVol',
                        u'turnoverValue',u'dealAmount', u'chgPct',
                        'turnoverRate',u'marketValue']):
    sql_str1 = sql_str_select_data1 % (','.join(field),ticker,begin,end)
    dataday = pd.read_sql(sql_str1,engine)
    dataday = dataday.applymap(lambda x: np.nan if x == 0 else x)
    dataday.rename(columns={'symbol':'ticker'},inplace=True)
    #升级后复权系数
    #后复权系数 MktEqudAdjAfGet
    sql_str_fq = """select tradeDate,accumAdjFactor from MktEqudAdjAfGet where 
    ticker = "%s" order by tradeDate"""
    sql_str_fq = sql_str_fq % ticker
    y = pd.read_sql(sql_str_fq,engine)
    dataday=pd.merge(dataday,y,on=['tradeDate'])
    return dataday.fillna(method = 'ffill')


## 得交易日历
def get_calender_range(begin, end):
    sql_str = """select tradeDate from yuqerdata.yq_index where symbol = "000001" 
    and tradeDate >="%s" and tradeDate <="%s" order by tradeDate""" % (begin, end)
    x=pd.read_sql(sql_str,engine)
    x=x['tradeDate'].values
    #b=[i.strftime('%Y-%m-%d') for i in x]
    return x

#获取所有交易日历
def get_calender():
    sql_str = '''select tradeDate from yuqerdata.yq_index where symbol = "000001" order by tradeDate'''
    x=pd.read_sql(sql_str,engine)
    x=x['tradeDate'].values
    #b=[i.strftime('%Y-%m-%d') for i in x]
    return x
#获取月度日历    
def get_month_calender(begin = '2000-01-01'):
    sql_str = '''select endDate from yuqerdata.yq_index_month where symbol = "000001" and endDate>="%s" order by endDate''' % (begin)
    x=pd.read_sql(sql_str,engine)
    x=x['endDate'].values
    #b=[i.strftime('%Y-%m-%d') for i in x]
    return x

#获取某个时间点的所有股票
def get_universe_date(tradeDate):
    sql_tmp = """select * from equget where listDate<'%s' and 
        (delistDate>'%s' or delistDate is null)"""
    x = pd.read_sql(sql_tmp % (tradeDate,tradeDate),engine)
    return x.ticker.tolist()
#获取A股所有的ticker
def get_symbol_A():
    sql_str = """select distinct(ticker)  from equget
                where equTypeCD = "A" and listStatusCD !="UN" and 
                ListSectorCD<=3 and length(ticker)=6  order by ticker"""
    x = pd.read_sql(sql_str,engine)
    return x.ticker.tolist()
#后复权月度数据
def get_MktEqumAdjGet(ticker,beginDate='2000-01-01', endDate='2049-01-01', field=u"ticker,endDate,chgPct"):      
    sql_str = """select %s from MktEqumAdjAfGet where endDate>="%s" and endDate<="%s" 
                order by endDate"""      
    sql_str = sql_str % (field,beginDate,endDate)
    #print(sql_str)
    return pd.read_sql(sql_str,engine)
#后复权月度收益 升级 截面同时选
def get_MktEqumAdjGet_update(beginDate='2000-01-01', endDate='2049-01-01', field=u"ticker,endDate,chgPct"):      
    sql_str = """select %s from MktEqumAdjAfGet where endDate>="%s" and endDate<="%s" 
                order by endDate"""      
    sql_str = sql_str % (field,beginDate,endDate)
    #print(sql_str)
    return pd.read_sql(sql_str,engine)
#get factor date
def get_MktStockFactorsOneDayGet(tradeDate,field=u""):
    if isinstance(field,list):
        filed_str = ','.join(field)
    else:
        filed_str = field
    sql_str = """select %s from yq_MktStockFactorsOneDayGet where tradeDate = "%s" order by tradeDate"""
    sql_str = sql_str % (filed_str,tradeDate)   
    #print(sql_str)
    return pd.read_sql(sql_str,engine)


def get_MktStockFactorsOneDayGet_special(factor_list,week_end_list):    
    sql_str_factor = 'select %s from yq_MktStockFactorsOneDayGet where tradeDate in (%s) '
    var_str = ','.join(['secID','ticker','tradeDate']+factor_list)
    date_str = '''"'''+'","'.join(week_end_list)+'''"'''
    sql_str_factor = sql_str_factor % (var_str,date_str)
    factor_df = pd.read_sql(sql_str_factor,engine)   
    return factor_df

#股票周行情
def get_MktEquwAdjAfGet(beginDate,endDate,field=r"*"):
    if isinstance(field,list):
        filed_str = ','.join(field)
    else:
        filed_str = field
    sql_str = """select %s from yq_MktEquwAdjAfGet where endDate >= "%s" 
                 and endDate<="%s" order by endDate """
    sql_str = sql_str % (filed_str,beginDate,endDate)   
    #print(sql_str)
    return pd.read_sql(sql_str,engine)    


def get_MktEqudAdjAfGet(beginDate,endDate,field=r""):
    if isinstance(field,list):
        filed_str = ','.join(field)
    else:
        filed_str = field
    sql_str = """select %s from yq_MktEqudAdjAfGet where tradeDate >= "%s" 
                 and tradeDate<="%s" order by tradeDate"""
    sql_str = sql_str % (filed_str,beginDate,endDate)   
    return pd.read_sql(sql_str,engine)

def get_MktEqudAdjAfGet_fill(beginDate,endDate,field=r""):
    if isinstance(field,list):
        filed_str = ','.join(field)
    else:
        filed_str = field
    sql_str = """select %s from MktEqudAdjAfGetF0S53 where tradeDate >= "%s" 
                 and tradeDate<="%s" order by tradeDate"""
    sql_str = sql_str % (filed_str,beginDate,endDate)   
    return pd.read_sql(sql_str,engine)

def get_MktEqudAdjAfGet_com(beginDate,endDate,field=r""):
    x0 = get_MktEqudAdjAfGet(beginDate,endDate,field)
    x1 = get_MktEqudAdjAfGet_fill(beginDate,endDate,field)
    x=pd.concat([x0,x1])
    x.reset_index(inplace=True,drop=True)
    return x
'''
def get_MktEqudAdjAfGet_com(beginDate,endDate,field=r""):
    if isinstance(field,list):
        filed_str = ','.join(field)
    else:
        filed_str = field
    sql_str = """select %s from MktEqudAdjAfGetF1S53 where tradeDate >= "%s" 
                 and tradeDate<="%s" order by tradeDate"""
    sql_str = sql_str % (filed_str,beginDate,endDate)   
    return pd.read_sql(sql_str,engine)
'''
    
def MktEqudGet(beginDate,endDate,field=r""):
    if isinstance(field,list):
        filed_str = ','.join(field)
    else:
        filed_str = field
    sql_str = """select %s from yq_dayprice where tradeDate >= "%s" 
                 and tradeDate<="%s" order by tradeDate"""
    sql_str = sql_str % (filed_str,beginDate,endDate)
    return pd.read_sql(sql_str,engine)
#filling=0部分
def MktEqudGet_fill(beginDate,endDate,field=r""):
    if isinstance(field,list):
        filed_str = ','.join(field)
    else:
        filed_str = field
    sql_str = """select %s from MktEqudGet0S53 where tradeDate >= "%s" 
                 and tradeDate<="%s" order by tradeDate"""
    sql_str = sql_str % (filed_str,beginDate,endDate)
    return pd.read_sql(sql_str,engine)

def MktEqudGet_com(beginDate,endDate,field=r""):
    x0=MktEqudGet(beginDate,endDate,field)
    x1=MktEqudGet_fill(beginDate,endDate,field)
    x =pd.concat([x0,x1])
    x.reset_index(inplace=True,drop=True)
    return x

    
def MktEqudGet_update(factor_list,week_end_list):
    if isinstance(factor_list,list):
        var_str = ','.join(factor_list)
    else:
        var_str = factor_list
    sql_str_factor = 'select %s from yq_dayprice where tradeDate in (%s) '    
    date_str = '''"'''+'","'.join(week_end_list)+'''"'''
    sql_str_factor = sql_str_factor % (var_str,date_str)
    x = pd.read_sql(sql_str_factor,engine)
    return x
def MktEqudGet_update_day(ticker_list,factor_list,begin,end):
    if isinstance(factor_list,list):
        var_str = ','.join(factor_list)
    else:
        var_str = factor_list
    sql_str_factor = 'select %s from yq_dayprice where symbol in (%s) and tradeDate>="%s" and tradeDate<="%s" '    
    ticker_list_str = '''"'''+'","'.join(ticker_list)+'''"'''
    sql_str_factor = sql_str_factor % (var_str,ticker_list_str,begin,end)
    x = pd.read_sql(sql_str_factor,engine)
    return x

    
def get_MktEqudAdjAfGet_update(ticker,begin,end,field=r""):
    if isinstance(field,list):
        filed_str = ','.join(field)
    else:
        filed_str = field
    sql_str = """select %s from yq_MktEqudAdjAfGet where ticker="%s" and tradeDate >= "%s" 
                 and tradeDate<="%s" order by tradeDate"""
    sql_str = sql_str % (filed_str,ticker,begin,end)   
    #print(sql_str)
    #print(sql_str)
    return pd.read_sql(sql_str,engine)

class time_use_tool():
    def __init__(self,ini_str=' '):
        self.t_now=time.time()
        self.t_now0=time.time()
        self.use('开始记录 %s ' % ini_str)
    def use(self,key_str=' ' ):
        tt = time.time()
        print('%s %s Time used %0.2f(All= %0.2f) seconds' % (key_str,time.ctime(),tt-self.t_now,tt-self.t_now0))
        self.t_now = tt

def get_week_month_tradeDate(start_date,end_date):
    #dates = pd.date_range(start_date,end_date,freq="D").astype(str)    
    # 获取月末交易日
    sql_str_calender = """select * from yuqerdata.yuqer_cal where exchangeCD = "%s" 
                        and calendarDate>="%s" and calendarDate<="%s" order by calendarDate"""    
    calendar_df = pd.read_sql(sql_str_calender % ('XSHG',start_date,end_date),engine)
    calendar_df['calendarDate'] = calendar_df['calendarDate'].astype(str)
    calendar_df['prevTradeDate'] = calendar_df['prevTradeDate'].astype(str)
    week_end_list = calendar_df[calendar_df['isWeekEnd']==1]['calendarDate'].values
    month_end_list = calendar_df[calendar_df['isMonthEnd']==1]['calendarDate'].values
    trade_date_list = calendar_df[calendar_df['isOpen']==1]['calendarDate'].values
    return trade_date_list,week_end_list,month_end_list

#时间函数
def get_week_month_tradeDate_update(start_date,end_date):
    #dates = pd.date_range(start_date,end_date,freq="D").astype(str)    
    # 获取月末交易日
    sql_str_calender = """select * from yuqerdata.yuqer_cal where exchangeCD = "%s" 
                        and calendarDate>="%s" and calendarDate<="%s" order by calendarDate"""    
    calendar_df = pd.read_sql(sql_str_calender % ('XSHG',start_date,end_date),engine)
    calendar_df['calendarDate'] = calendar_df['calendarDate'].astype(str)
    calendar_df['prevTradeDate'] = calendar_df['prevTradeDate'].astype(str)
    week_end_list = sorted(calendar_df[calendar_df['isWeekEnd']==1]['calendarDate'].tolist())
    month_end_list = sorted(calendar_df[calendar_df['isMonthEnd']==1]['calendarDate'].tolist())
    trade_date_list = sorted(calendar_df[calendar_df['isOpen']==1]['calendarDate'].tolist())
    #daily_trade_list = sorted(cal_dates_df.query("isOpen==1")['calendarDate'].tolist())
    return trade_date_list,week_end_list,month_end_list,calendar_df

#S49添加 根据endDate选择数据
def get_FdmtBSGet(beginDate,endDate,field):
    sql_str_factor = 'select %s from yq_FdmtBSGet where endDate >="%s" and endDate<="%s" '
    if isinstance(field,list):
        var_str = ','.join(field)
    else:
        var_str = field
    sql_str_factor = sql_str_factor % (var_str,beginDate,endDate)
    factor_df = pd.read_sql(sql_str_factor,engine)   
    return factor_df

def get_FdmtISQPITGet(beginDate,endDate,field):
    sql_str_factor = 'select %s from yq_FdmtISQPITGet where endDate >="%s" and endDate<="%s" '
    if isinstance(field,list):
        var_str = ','.join(field)
    else:
        var_str = field
    sql_str_factor = sql_str_factor % (var_str,beginDate,endDate)
    factor_df = pd.read_sql(sql_str_factor,engine)   
    return factor_df

def FdmtEeGet(beginDate,endDate,field,key_words = 'publishDate'):
    if isinstance(field,list):
        var_str = ','.join(field)
    else:
        var_str = field
    sql_str_factor = 'select %s from yq_FdmtEeGet where %s >="%s" and %s<="%s" '
    sql_str_factor = sql_str_factor % (var_str,key_words,beginDate,key_words,endDate)
    factor_df = pd.read_sql(sql_str_factor,engine)   
    return factor_df

def FdmtEfGet(beginDate,endDate,field,key_words = 'publishDate'):
    if isinstance(field,list):
        var_str = ','.join(field)
    else:
        var_str = field
    sql_str_factor = 'select %s from yq_FdmtEfGet where %s >="%s" and %s<="%s" '
    sql_str_factor = sql_str_factor % (var_str,key_words,beginDate,key_words,endDate)
    factor_df = pd.read_sql(sql_str_factor,engine)   
    return factor_df
#FdmtISGet
def FdmtISGet(beginDate,endDate,field,key_words = 'publishDate'):
    if isinstance(field,list):
        var_str = ','.join(field)
    else:
        var_str = field
    sql_str_factor = 'select %s from nincome where %s >="%s" and %s<="%s" '
    sql_str_factor = sql_str_factor % (var_str,key_words,beginDate,key_words,endDate)
    factor_df = pd.read_sql(sql_str_factor,engine)   
    return factor_df    

def get_a_stock_pool():
    sql_str_temp = """select ticker from equget where equTypeCD="A" and listStatusCD in
    ("L","S","DE")"""
    x=pd.read_sql(sql_str_temp,engine)['ticker'].tolist()
    x.remove('DY600018')
    return x

#获取行业分类
def get_industry_class(t):
    sql_str1 = '''select ticker,industryID1 from yuqerdata.yq_industry where 
                industryVersionCD="010303" and intodate <= "%s" and 
                (outDate>"%s" or outDate is null)''' % (t,t)
    x = pd.read_sql(sql_str1,engine)
    return x
#获取行业分类
def get_industry_class_update(t,fileds_str = 'ticker,industryID1'):
    sql_str1 = '''select %s from yuqerdata.yq_industry where 
                industryVersionCD="010303" and intodate <= "%s" and 
                (outDate>"%s" or outDate is null)''' % (fileds_str,t,t)
    x = pd.read_sql(sql_str1,engine)
    return x
#
def get_file_name(file_dir,file_type):
    L=[]
    L_s = []   
    for root, dirs, files in os.walk(file_dir):  
        for file in files:  
            if os.path.splitext(file)[1] == file_type:  
                L.append(os.path.join(root, file))  
                L_s.append(file)
    return L,L_s

def get_ini_data(tn,var_name,db):
    sql_str = 'select %s from %s order by %s desc limit 1' % (var_name,tn,var_name)
    t0 = pd.read_sql(sql_str,db)
    return t0[var_name].astype(str).values[0]

#两个list并集
def list_com(a,b):
    return list(set(a).union(set(b)))

#两个list交集
def list_intersec(a,b):
    return list(set(a).intersection(set(b)))

#两个list差
def list_diff(a,b):
    return list(set(a).difference(set(b)))

def save_pickle(fn,x):
    with open(fn, 'wb') as f:
        pickle.dump(x, f)
def read_pickle(fn):    
    with open(fn, 'rb') as f:
        return pickle.load(f)
#载入yuqer因子
def get_factor_data_update(fac_name,target_date,method='outer'):
    def fac_info_trans(id_sec):
        factor_style10_str0 = 'select symbol,f_val as %s from %s where tradingdate ="%s"'
        factor_style10_str1 = 'select ticker as symbol,%s from %s where tradeDate ="%s"'
        if id_sec==0:
            sub_tn = db_name_factor0
            factor_style10_str = factor_style10_str1
            sub_engine = engine
        elif id_sec==-1:
            sub_tn = db_name_factor1
            factor_style10_str = factor_style10_str1
            sub_engine = engine
        else:
            sub_tn = ''
            factor_style10_str = factor_style10_str0
            sub_engine = engine_yq_cub
        return sub_tn,factor_style10_str,sub_engine
    #def get_factor_data_update(fac_name,target_date):
    t_cut = '2020-07-01'
    db_name_factor0 = 'yq_MktStockFactorsOneDayGet'
    db_name_factor1 = 'yq_MktStockFactorsOneDayProGet'
    fac_name = [i.lower() for i in fac_name]
    
    tn = ['yq_MktStockFactorsOneDayGet','yq_MktStockFactorsOneDayProGet']
    var = [pd.read_sql('desc %s ' % i,engine).Field.tolist() for i in tn]
    var[0] = [i.lower() for i in var[0]]
    var[1] = [i.lower() for i in var[1]]
        
    fac_id = [0 if i in var[0] else 1 if i in var[1] else 3 for i in fac_name]
    factor_style10=fac_name.copy()
    factor_style_tn = fac_id.copy()
    x=pd.DataFrame()
    
    obj_time = time_use_tool('载入数据立方因子 %s' % target_date)
    fac0=[]
    fac1=[]
    for i,j in enumerate(factor_style_tn):
        if j==0:
            fac0.append(factor_style10[i])
        elif j==1:
            fac1.append(factor_style10[i])
    
    if len(fac0)>0:
        sub_tn,factor_style10_str,sub_engine=fac_info_trans(0)
        sub_sql_str = factor_style10_str % (','.join(fac0),sub_tn,target_date)
        sub_x0 = pd.read_sql(sub_sql_str,sub_engine)
    else:
        sub_x0= pd.DataFrame()
    obj_time.use('loading factor %d-%d out of %d' % (0,len(fac0),len(factor_style_tn)))   
    sub_x=[]
    sub_y=[]
    if len(fac1)>0:
        if target_date>=t_cut:
            sub_tn,factor_style10_str,sub_engine=fac_info_trans(-1)
            sub_sql_str = factor_style10_str % (','.join(fac1),sub_tn,target_date)
            sub_y = pd.read_sql(sub_sql_str,sub_engine)
        else:        
            T = len(fac1)
            _,factor_style10_str,sub_engine=fac_info_trans(1)
            sub_y=pd.DataFrame()
            for i,sub_f_name in enumerate(fac1):
                sub_tn = sub_f_name.lower()
                sub_sql_str = factor_style10_str % (sub_f_name,sub_tn,target_date)
                sub_x = pd.read_sql(sub_sql_str,sub_engine)
                if len(sub_x)==0:
                    continue
                if len(sub_y)==0:
                    sub_y = sub_x
                else:
                    sub_y = pd.merge(sub_y,sub_x,on='symbol',how = method)
                obj_time.use('loading factor %s %d out of %d' % (sub_f_name,i,T))
        
    if len(sub_x0)>0 and len(sub_y)==0:
        x=sub_x0
    elif len(sub_x0)==0 and len(sub_y)>0:
        x=sub_y
    else:
        x=pd.merge(sub_x0,sub_y,on='symbol',how = method)
    return x

def get_factor_data_update2(fac_name,target_date,method='outer'):
    factor_df=get_factor_data_update(fac_name,target_date,method)
    c={}
    for sub_f in fac_name:
        c[sub_f.lower()]=sub_f
    factor_df.rename(columns=c,inplace=True)
    factor_df.rename(columns={'symbol':'ticker'},inplace=True)
    factor_df['secID'] = factor_df.ticker.apply(lambda x:ticker2secID(x))
    return factor_df
    
def get_factor_update(field,week_end_list):
    def get_factor_thread(inputdata):
        x,y=inputdata
        f = get_factor_data_update(x,y)
        f['tradeDate'] = y
        return f
    pool = ThreadPool(processes=num_core)
    fac_name=[field for i in week_end_list]
    temp = pool.map(get_factor_thread, zip(fac_name,week_end_list))
    pool.close()
    pool.join() 
    factor_df=pd.concat(temp)
    #名称对齐
    c={}
    for sub_f in field:
        c[sub_f.lower()]=sub_f
    factor_df.rename(columns=c,inplace=True)
    factor_df.rename(columns={'symbol':'ticker'},inplace=True)
    factor_df['secID'] = factor_df.ticker.apply(lambda x:ticker2secID(x))
    return factor_df

def add_0(x):
    if isinstance(x,int):
        x= '%0.6d' % x
    else:
        x=x.rjust(6,'0')
    return x
#ST info
def SecSTGet(t0,t1):
    sql_tmp = 'select * from yuqerdata.st_info where tradeDate>="%s" and tradeDate<="%s"'
    ticker=  pd.read_sql(sql_tmp % (t0,t1),engine)
    ticker['ticker'] = ticker.ticker.apply(lambda x:add_0(x))
    return ticker

def ticker2secID(ticker):
    """
    ticker转换secID
    转换规则：secID = ticker + 后缀：如果股票属于沪市，则后缀为'.XSHG'，如果属于深市，则后缀为'.XSHE'
    """
    ticker = '0'*(6-len(ticker)) + ticker
    if ticker[0] == '6':
        secID = ticker + '.XSHG'
    else:
        secID = ticker + '.XSHE'
    return secID

def get_trade_dates(start_date, end_date, frequency='d'):
    """
    输入起始日期和频率，即可获得日期列表（daily包括起始日，其余的都是位于起始日中间的）
    输入：
       start_date，开始日期，'YYYYMMDD'形式
       end_date，截止日期，'YYYYMMDD'形式
       frequency，频率，daily为所有交易日，weekly为每周最后一个交易日，monthly为每月最后一个交易日，quarterly为每季最后一个交易日
    返回：
       获得list型日期列表，以'YYYYMMDD'形式存储
    """
    #data = DataAPI.TradeCalGet(exchangeCD=u"XSHG", beginDate=start_date, endDate=end_date,
    #                           field=u"calendarDate,isOpen,isWeekEnd,isMonthEnd,isQuarterEnd", pandas="1")
    field=u"calendarDate,isOpen,isWeekEnd,isMonthEnd,isQuarterEnd"
    #yuqer_cal
    sql_tmp = """select %s from yuqer_cal where exchangeCD = "XSHG" and calendarDate>="%s" 
                and calendarDate<= "%s" order by calendarDate""" % (field,start_date,end_date)
    data = pd.read_sql(sql_tmp,engine)
    if frequency == 'd':
        data = data[data['isOpen'] == 1]
    elif frequency == 'w':
        data = data[data['isWeekEnd'] == 1]
    elif frequency == 'm':
        data = data[data['isMonthEnd'] == 1]
    elif frequency == 'q':
        data = data[data['isQuarterEnd'] == 1]
    else:
        raise ValueError('调仓频率必须为d/w/m！！')
    #date_list = map(lambda x: x[0:4] + x[5:7] + x[8:10], data['calendarDate'].values.tolist())
    date_list = data['calendarDate'].astype(str).values.tolist()
    return date_list

def month_tool(x,N,method='today'):
    if isinstance(x,str):
        x = datetime.strptime(x,'%Y-%m-%d')
    year = x.year
    mon  = x.month
    day  = x.day
    hour = x.hour
    minu  = x.minute
    sec  = x.second
    def today():
        '''
        get today,date format="YYYY-MM-DD"
        '''
        return date.today()
    def todaystr():
        '''
        get date string
        date format="YYYYMMDD"
        '''
        return year+mon+day
    #def datetime():
    #    '''
    #    get datetime,format="YYYY-MM-DD HH:MM:SS"
    #    '''
    #    return strftime("%Y-%m-%d %H:%M:%S",localtime())
    def datetimestr():
        '''
        get datetime string
        date format="YYYYMMDDHHMMSS"
        '''
        return year+mon+day+hour+minu+sec
    def getdayofday(n=0):
        '''
        if n>=0,date is larger than today
        if n<0,date is less than today
        date format = "YYYY-MM-DD"
        '''
        if(n<0):
            n = abs(n)
            return date.today()-timedelta(days=n)
        else:
            return date.today()+timedelta(days=n)
    def getdaysofmonth(year,mon):
        '''
        get days of month
        '''
        return calendar.monthrange(year, mon)[1]
    def getfirstdayofmonth(year,mon):
        '''
        get the first day of month
        date format = "YYYY-MM-DD"
        '''
        days="01"
        if(int(mon)<10):
            mon = "0"+str(int(mon))
        arr = (year,mon,days)
        return "-".join("%s" %i for i in arr)
    def getlastdayofmonth(year,mon):
        '''
        get the last day of month
        date format = "YYYY-MM-DD"
        '''
        days=calendar.monthrange(year, mon)[1]
        mon = addzero(mon)
        arr = (year,mon,days)
        return "-".join("%s" %i for i in arr)
    def get_firstday_month(n=0):
        '''
        get the first day of month from today
        n is how many months
        '''
        (y,m,d) = getyearandmonth(n)
        d = "01"
        arr = (y,m,d)
        return "-".join("%s" %i for i in arr)
    def get_lastday_month(n=0):
        '''
        get the last day of month from today
        n is how many months
        '''
        return "-".join("%s" %i for i in getyearandmonth(n))
        
    def get_today_month(n=0):
        '''
        get last or next month's today
        n is how many months
        date format = "YYYY-MM-DD"
        '''
        (y,m,d) = getyearandmonth(n)
        arr=(y,m,d)
        if(int(day)<int(d)):
            arr = (y,m,day)
        return "-".join("%s" %i for i in arr)
    def getyearandmonth(n=0):
        '''
        get the year,month,days from today
        befor or after n months
        '''
        thisyear = int(year)
        thismon = int(mon)
        totalmon = thismon+n
        if(n>=0):
            if(totalmon<=12):
                days = str(getdaysofmonth(thisyear,totalmon))
                totalmon = addzero(totalmon)
                return (year,totalmon,days)
            else:
                i = math.floor(totalmon/12)
                j = totalmon-i*12
                thisyear += i
                days = str(getdaysofmonth(thisyear,j))
                j = addzero(j)
                return (str(thisyear),str(j),days)
        else:
            if((totalmon>0) and (totalmon<12)):
                days = str(getdaysofmonth(thisyear,totalmon))
                totalmon = addzero(totalmon)
                return (year,totalmon,days)
            else:
                i = -math.ceil(abs(totalmon)/12)
                j = totalmon -i*12
                if(j==0):
                    i-=1
                    j=12
                thisyear +=i
                days = str(getdaysofmonth(thisyear,j))
                j = addzero(j)
                return (str(thisyear),str(j),days)
    def addzero(n):
        '''
        add 0 before 0-9
        return 01-09
        '''
        nabs = abs(int(n))
        if(nabs<10):
            return "0"+str(nabs)
        else:
            return nabs
    if method == 'today' or method =='t':
        return get_today_month(N)
    elif method == 'lastday' or method=='l':
        return get_lastday_month(N)
    else:
        return get_firstday_month(N)
#回退数据日期
def get_back_date(date,N,trade_list):
    date_id = trade_list.index(date)-N
    if date_id<0:
        date_id = 0
    elif date_id>=len(trade_list)-1:
        date_id = len(trade_list)-1
    return trade_list[date_id]
#成分股
def get_IdxConsGet(ticker,intoDate):
    if len(ticker)==6:
        sql_str1 = 'select tradeDate from IdxConsGetS57 where tradeDate<="%s" and ticker="%s" order by tradeDate desc limit 1'
        t0 = pd.read_sql(sql_str1 % (intoDate,ticker),engine)
        t0 = t0.tradeDate.astype(str).tolist()[0]
        sql_strl = 'select * from IdxConsGetS57 where ticker="%s" and tradeDate="%s"'
        x=pd.read_sql(sql_strl % (ticker,t0),engine)
    else:
        sql_str1 = 'select tradeDate from IdxConsGetS57 where tradeDate<="%s"  order by tradeDate desc limit 1'
        t0 = pd.read_sql(sql_str1 % (intoDate),engine)
        t0 = t0.tradeDate.astype(str).tolist()[0]
        sql_strl = 'select * from IdxConsGetS57 where tradeDate="%s"'
        x=pd.read_sql(sql_strl % (t0),engine)
    return x.consTickerSymbol.unique(),x.consID.unique()
    