# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 20:29:30 2021
更新最新月份数据，至少需要前12个月的数据才可以
@author: ASUS
"""
# coding=utf-8

from tqdm import tqdm
import pandas as pd
import numpy as np
import math
from yq_toolsS45 import get_week_month_tradeDate_update
import quant_utilS73 as qutil
import datetime
from sqlalchemy.types import NVARCHAR, Float,DATE,Integer
from yq_toolsS45 import engine
from yq_toolsS45 import MktIdxmGet
from yq_toolsS45 import time_use_tool
from yq_toolsS45 import get_MktEqumAdjGet_update
from yq_toolsS45 import get_industry_class_update as EquIndustryGet
from yq_toolsS45 import RMExposureDayGet
from yq_toolsS45 import get_factor_update
from yq_toolsS45 import create_db
from yq_toolsS45 import get_table_date
from yq_toolsS45 import table_in_database

tn_f1 = 'symbol_pool_s76p1'
eg37 = create_db('s37')


obj_t = time_use_tool()
rt_key1 = 'S76Part1-'

def report(info,rt_key = rt_key1):
    obj_t.use('%s-%s' % (rt_key,info))

def trans_date(x):
    if len(x)==8:
        x = '%s-%s-%s' % (x[:4],x[4:6],x[6:])
    else:
        x = x.replace('-','')
    return x


# 基础函数准备
def fac_process(factor_df, factor_list):
    '''
    因子处理函数
    '''
    # 去极值
    w_factor_df = qutil.mad_winsorize(factor_df, factor_list, sigma_n=3)
    
    # 完全中性化
    n_factor_df = qutil.neutralize_dframeV2(w_factor_df.copy(), factor_list, exclude_style=[])

    # 标准化
    s_factor_df = n_factor_df.copy()
    s_factor_df[factor_list] = s_factor_df.groupby('tradeDate')[factor_list].apply(lambda df: (df-df.mean()) / df.std())
    return s_factor_df



def annual_std_fac(df, fac_col_name, std_col_name, span=12):
    '''
    函数：相对因子过去一年数据进行标准化处理，求解动态因子
    参数：
        df: DataFrame，某只股票时间序列因子值，columns至少包含tradeDate、fac_col_name
        fac_col_name: 待标准化处理因子在df中列名
        std_col_name: 标准化处理后的因子列名    
        span: 滚动历史数据的期数
    返回：DataFrame，columns至少包含tradeDate、fac_col_name、std_col_name
    '''
    df['mean'] = df[fac_col_name].rolling(window=span, min_periods=10).mean() 
    df['std'] = df[fac_col_name].rolling(window=span, min_periods=10).std()
    df['std'] = np.where(df['std']==0, np.nan, df['std'])
    df[std_col_name] = (df[fac_col_name]-df['mean']) / df['std']
    return df[[std_col_name]]


if __name__ == "__main__":
    report('开始')
    #初始时间
    if table_in_database('s37',tn_f1):
        t0=get_table_date(tn_f1,eg37,'tradeDate')
        sdate = datetime.datetime.strptime(t0,'%Y-%m-%d')-datetime.timedelta(days=(366+61))
        sdate = sdate.strftime('%Y-%m-%d')    
    else:
        t0 = '2009-01-01'
        sdate = t0
    
    edate = '2099-01-01'
    
    report('载入基本数据')
    trade_date_list,week_end_list,month_end_list,cal_dates_df = get_week_month_tradeDate_update(sdate,edate)
    trade_date_list = [trans_date(i) for i in trade_date_list]
    week_end_list = [trans_date(i) for i in week_end_list]
    month_end_list = [trans_date(i) for i in month_end_list] 
    if max(month_end_list)>t0.replace('-',''):
        cal_dates_df.calendarDate = cal_dates_df.calendarDate.astype(str).apply(lambda x:trans_date(x))
        cal_dates_df.prevTradeDate = cal_dates_df.prevTradeDate.astype(str).apply(lambda x:trans_date(x))    
        # 全A投资域
        sql_str_temp = """select ticker from equget where equTypeCD="A" and listStatusCD in
            ("L")"""
        stock_list=pd.read_sql(sql_str_temp,engine)['ticker'].tolist()
        stock_list = [i for i in stock_list if len(i)==6]
            
        forbidden_pool = qutil.stock_special_tag(sdate, edate, pre_new_length=60)  # 次新股、st股、停牌个股
        # 筛选一字板个股
        #手动筛选
        sql_tmp = 'select tradeDate,symbol as ticker from yq_dayprice where highestPrice=lowestPrice and highestPrice>0 and tradeDate>="%s" and tradeDate<="%s"'
        mkt_df = pd.read_sql(sql_tmp % (sdate,edate),engine)
        mkt_df = mkt_df[mkt_df.ticker.isin(stock_list)]
        mkt_df['tradeDate'] = mkt_df['tradeDate'].astype(str).apply(lambda x: x.replace('-', ''))
        limit_df = mkt_df.copy()
        limit_df['special_flag'] = 'limit'
        forbidden_pool = forbidden_pool.append(limit_df)
        forbidden_pool = forbidden_pool.merge(cal_dates_df, left_on=['tradeDate'], right_on=['calendarDate'])
        forbidden_pool = forbidden_pool[['ticker', 'tradeDate', 'prevTradeDate', 'special_flag']]
    
    
        report("风格因子计算...")
        style_list = ['BETA', 'SIZE', 'RESVOL', 'MOMENTUM', 'SIZENL', 'EARNYILD', 'BTOP', 'GROWTH', 'LEVERAGE', 'LIQUIDTY']
        field=['ticker','tradeDate'] + style_list
        style_fac_df = RMExposureDayGet(sdate,edate,field)
        style_fac_df.tradeDate = style_fac_df.tradeDate.astype(str).apply(lambda x:x.replace('-',''))
        style_fac_df.sort_values(['ticker', 'tradeDate'], inplace=True)
        
        report("市盈率、历史归母净利润增速和流通市值数据")    
        key_str = ['PE', 'NPParentCompanyGrowRate', 'NegMktValue']
        pe_df = get_factor_update(key_str,month_end_list)
        pe_df = pe_df.sort_values(['ticker', 'tradeDate'])
        
        report("一致预期净利润数据")       
        field="secCode,repForeTime,ConProfitType,foreYear,conProfit"
        sql_tmp = 'select %s from ResConSecDataGet18 where repForeTime = "%s"'
        con_profit_list = []
        for td in tqdm(month_end_list):
            tmp_df = pd.read_sql(sql_tmp % (field,td),engine)
            con_profit_list.append(tmp_df)
        con_profit_df = pd.concat(con_profit_list)
        con_profit_df['repForeTime'] = con_profit_df['repForeTime'].astype(str).apply(lambda x: x[:10].replace('-', ''))
        con_profit_df = con_profit_df.sort_values(['secCode', 'repForeTime', 'foreYear'])
        # 历史净利润增速计算BET和PEG因子
        pe_df['PEG_h'] = np.where(pe_df.NPParentCompanyGrowRate==0, np.nan, pe_df['PE']/pe_df['NPParentCompanyGrowRate'])
        pe_df['PEG_h'] = np.where((pe_df.PE<0)|(pe_df.NPParentCompanyGrowRate<0), np.nan, pe_df['PEG_h']) 
        pe_df['BET_h'] = np.log((1+pe_df['NPParentCompanyGrowRate']*(pe_df['PE'])).tolist()) / np.log((1+pe_df['NPParentCompanyGrowRate']).tolist())
        pe_df['BET_h'] = np.where((pe_df['PE']<0)|(pe_df['NPParentCompanyGrowRate']<=0), np.nan, pe_df['BET_h']) 
    
        report("企业成长性估计：计算不同类型的归母净利润增长率：")
        # 剔除ConProfitType=3的字段
        con_profit_df = con_profit_df.query("ConProfitType<3").sort_values(['secCode', 'repForeTime', 'foreYear']) 
        # 提取基期数据
        baseyr_df = con_profit_df.query("ConProfitType==0").drop_duplicates(['secCode', 'repForeTime'], keep='last')
        # 提取一致预期数据
        forecast_df = con_profit_df.query("ConProfitType>0").reset_index(drop=True)
        # 最近一期的一致预期净利润数据
        min_forecast_df = forecast_df.drop_duplicates(['secCode', 'repForeTime'], keep='first')
        # 最近第二期的一致预期净利润数据
        second_forecast_df = forecast_df.drop(min_forecast_df.index, axis=0).drop_duplicates(['secCode', 'repForeTime'], keep='first')
        
        merge_df = baseyr_df.merge(min_forecast_df, on=['secCode', 'repForeTime'], how='left', suffixes=('', '_min'))
        merge_df = merge_df.merge(second_forecast_df, on=['secCode', 'repForeTime'], how='left', suffixes=('', '_sec'))
        merge_df['FY1_gap'] = np.where(merge_df['conProfit_min']!=0, merge_df['foreYear_min'] - merge_df['foreYear'], np.nan)
        merge_df['FY2_gap'] = np.where(merge_df['conProfit_sec']!=0, merge_df['foreYear_sec'] - merge_df['foreYear'], np.nan)
        # 计算一致预期FY2相对历史真实值的增长率（预测期超过1期）
        merge_df['cagr_growth'] = merge_df.apply(lambda s: np.sign(s['conProfit_sec']) * \
                math.pow(abs(s['conProfit_sec']/s['conProfit']), (1.0/2))-1 if s['FY2_gap']==2 else np.nan, axis=1)
        # 计算一致预期FY1相对历史真实值的增长率（预测期仅1期）
        merge_df['fy1_growth'] = np.where(merge_df['FY1_gap']==1, (merge_df['conProfit_min']-merge_df['conProfit']) / abs(merge_df['conProfit']), np.nan)
        fac_df = merge_df[['secCode', 'repForeTime', 'cagr_growth', 'fy1_growth', 'FY2_gap', 'FY1_gap'\
                           ]].rename(columns={'secCode':'ticker', 'repForeTime':'tradeDate'}).merge(pe_df, on=['ticker', 'tradeDate'])
           
        # 计算不同类型的BET和PEG因子
        fac_df['BET_c'] = np.log((1+fac_df['cagr_growth']*(fac_df['PE'])).tolist()) / np.log((1+fac_df['cagr_growth']).tolist())  
        fac_df['BET_c'] = np.where((fac_df['PE']<0)|(fac_df['cagr_growth']<=0), np.nan, fac_df['BET_c'])
        fac_df['PEG_c'] = np.where(fac_df.cagr_growth==0, np.nan, fac_df['PE']/fac_df['cagr_growth'])
        fac_df['PEG_c'] = np.where((fac_df.PE<0)|(fac_df.cagr_growth<0), np.nan, fac_df['PEG_c']) 
    
        report("计算时间序列标准化的动态因子：")
        fac_df.sort_values(['ticker','tradeDate'], inplace=True)
        fac_df['PEG_c_STD'] = fac_df.groupby('ticker').apply(lambda df: annual_std_fac(df, 'PEG_c', 'PEG_c_STD'))
        fac_df['BET_c_STD'] = fac_df.groupby('ticker').apply(lambda df: annual_std_fac(df, 'BET_c', 'BET_c_STD'))
        factor_list = ['BET_c', 'BET_c_STD', 'PEG_c', 'PEG_c_STD']
        
        report('因子中性化')
        s_fac_df = fac_df.copy()
        s_fac_df[factor_list] = s_fac_df[factor_list].apply(lambda x: x.astype(float))
        s_fac_df[factor_list] = - s_fac_df[factor_list]  # 因子方向调整
        
        # 完全中性化 太耗时
        a_factor_df = fac_process(s_fac_df, factor_list)
        a_factor_df = a_factor_df[['ticker','tradeDate']+factor_list]
        a_factor_df.dropna(subset = ['BET_c','BET_c_STD','PEG_c','PEG_c_STD'],how='all',inplace=True)
        a_factor_df = a_factor_df[a_factor_df.tradeDate>t0.replace('-','')]
        var = ['ticker', 'tradeDate', 'BET_c', 'BET_c_STD', 'PEG_c', 'PEG_c_STD']
        var_t = [NVARCHAR(6),DATE,Float,Float,Float,Float]
        dtypedict1 = dict(zip(var,var_t))
        a_factor_df.to_sql(tn_f1,eg37,if_exists='append',index=False,chunksize=3000,dtype=dtypedict1)
        report('完成')
    else:
        report('因子数据已经是最新，无需更新')