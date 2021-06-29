#utf-8
"""
精简策略
pandas 1.1.4
对接S37程序
20210516
扩展计算所有A股的结果，便于后续分析
"""

import numpy as np
from alg_macd_direction import macd_direction
import pandas as pd
from yq_toolsS45 import engine
from yq_toolsS45 import save_pickle
from yq_toolsS45 import get_MktEqudAdjAfGet_update
from yq_toolsS45 import table_in_database
from yq_toolsS45 import create_db
from yq_toolsS45 import MktIdxdGet
from sqlalchemy.types import NVARCHAR, Float,DATE,Integer
from yq_toolsS45 import time_use_tool
from yq_toolsS45 import create_table
from yq_toolsS45 import get_symbol_A
from yq_toolsS45 import get_MktEqudAdjAfGet
from yq_toolsS45 import get_table_date
from yq_toolsS45 import get_spx_com
import datetime

from tqdm import tqdm
import multiprocessing
num_core2 = int(multiprocessing.cpu_count()/2)

BX_sel = True  #False#并行选择
disp = False
method_sel = 1
obj_t = time_use_tool()

dn = 'parapool'
eg = create_db(dn)
eg_pro = create_db('data_pro','localhost')
engine37 = create_db('s37','localhost')
if method_sel==2:
    tn = 's76p2_signal'
else:
    tn = 's76p2_signal_signal1'
if disp:
    import matplotlib.pyplot as plt
tn_var = ['ticker', 'dtype','tradeDate', 'CHGPct', 'sig']
t00 = '2010-01-01'
dtypedict1 = dict(zip(tn_var,[NVARCHAR(10),NVARCHAR(10),DATE,Float,Integer]))

fee1 = dict(zip(['US','HK','forex_day','as51','topix','twse','csi_stock','hsce','hk_ggt'],[1,11,0.25,1,1,2,2,11,11]))
fee2 = dict(zip(['US','HK','forex_day','as51','topix','twse','csi_stock','hsce','hk_ggt'],[1,11,0.25,1,1,32,12,11,11]))




index_id_zx02 = ['kosdaq', 'kospi', 'msci', 'ndx', 'nifty', 'nky', 'RTY', 'set50', 'sx5e',
                       'ukx', 'xin9i']
index_id_zx02_info = ['KOSDAQ','KOSPI2','TAMSCI','NDX','NIFTY','NKY','RTY','SET50','SX5E',
                      'UKX','XIN9I']
fee1_tmp = [3/10000,3/10000,2/10000,1/10000,1/1000,1/10000,1/10000,11/10000,1/10000,1/10000,2/10000]
fee2_tmp = [3/1000,3/1000,32/10000,1/10000,1/1000,1/10000,1/10000,22/10000,1/10000,1/10000,12/10000]
index_tdx = {'as51':'AS51','topix':'TPX','twse':'TWSE','hsce':'HSCEI','hk_ggt':'hk_ggt'}
index_tdx.update(dict(zip(index_id_zx02,index_id_zx02_info)))
tn_tdx = 'main_index_s68'
#tn_tdx = 'main_index_zx02'
#手续费
fee1_1 = dict(zip(index_id_zx02,[i*100 for i in fee1_tmp]))
fee2_1 = dict(zip(index_id_zx02,[i*100 for i in fee2_tmp]))
fee1.update(fee1_1)
fee2.update(fee2_1)



def get_ticker_pool(dtype):
    if dtype=="csi":
        hs_300_pool = get_symbol_A()
    elif dtype == "US":
        hs_300_pool = get_spx_com('2021-05-09')    
    elif dtype == "HK":
        hs_300_pool = ' select distinct(ticker) from main_index_s68 where index_id = "HSI" and ticker != "HSI"'
        hs_300_pool = pd.read_sql(hs_300_pool,eg_pro)
        hs_300_pool = hs_300_pool.ticker.tolist()
    elif dtype == 'hsce':
        hs_300_pool = ' select distinct(ticker) from main_index_s68 where index_id = "hsce" and ticker != "HSCEI"'
        hs_300_pool = pd.read_sql(hs_300_pool,eg_pro)
        hs_300_pool = hs_300_pool.ticker.tolist()
        hs_300_pool = ['%0.5d' % int(i) for i in hs_300_pool]
    elif dtype == 'hk_ggt':
        hs_300_pool = pd.read_csv('港股通股票池20210513.csv',encoding='gbk')
        hs_300_pool = hs_300_pool.ticker.tolist()
        hs_300_pool = ['%0.5d' % i for i in hs_300_pool]
        
        tmp1 = get_ticker_pool('HK')
        tmp2 = get_ticker_pool('hsce')
        hs_300_pool = list(set(hs_300_pool+tmp1+tmp2))
    elif dtype in index_tdx.keys():
        hs_300_pool = ' select distinct(ticker) from %s where index_id = "%s" and ticker != "%s"'
        hs_300_pool = pd.read_sql(hs_300_pool % (tn_tdx,dtype,index_tdx[dtype]),eg_pro)
        hs_300_pool = hs_300_pool.ticker.tolist()    
    elif dtype in ["forex_day","forex_min"]:
        hs_300_pool = ['EURUSD', 'GBPUSD', 'USDJPY', 'AUDUSD', 'NZDUSD', 'USDCAD',
               'USDCHF', 'EURGBP', 'EURCHF', 'EURAUD', 'EURNZD', 'EURCAD', 'AUDJPY', 'NZDJPY',
               'EURJPY', 'CHFJPY', 'GBPJPY', 'CADJPY', 'USDCNH', 'USDINR', 'USDTRY', 'USDRUB',
               'USDZAR'] 
    else:
        hs_300_pool = []
    return hs_300_pool
def ini_table_76p2(tn):
    var_name =  ['ticker', 'dtype','tradeDate', 'CHGPct', 'sig','sig1','fee']
    var_type = []
    for sub_info in var_name:
        if sub_info in ["ticker",'dtype']:
            var_type.append('varchar(10)')
        elif sub_info == 'tradeDate':
            var_type.append('date')
        elif sub_info in ['sig','sig1']:
            var_type.append('int')
        else:
            var_type.append('float')
    key_str = 'ticker,tradeDate'            
    create_table(dn,tn,var_name,var_type,key_str)  

def get_signal(s,macd_dir):
    break_flag = True
    # macd的区间划分不存在未来信息
    while break_flag:
        #print('ok', macd_dir.shape[0], s.shape[0])
        max_dir, min_dir = np.zeros_like(macd_dir), np.zeros_like(macd_dir)
        max_last, min_last = np.zeros_like(macd_dir), np.zeros_like(macd_dir)
        # init
        confirm_up, confirm_down, exp_flag = False, False, False
        b_flag, c_flag = False, False
        b_ind = 0
        t_max, t_min, t_dir = 1, 10e6, macd_dir[0]
        l_max, l_min = s[0], s[0]
        for ind, (s_, dir_) in enumerate(zip(s.values, macd_dir)):
            if t_dir != dir_:
                l_max, l_min = t_max, t_min
                t_max, t_min = s_, s_
                l_ind = b_ind
                b_ind = ind
    
            # update min max
            if t_max < s_:
                t_max = s_
                # update before max info
                max_dir[b_ind:ind + 1] = t_max
            if t_min > s_:
                t_min = s_
                # update before min info
                min_dir[b_ind:ind + 1] = t_min
    
            # update exp
            if dir_ == 1 and t_min < l_min and t_max < l_max:  # up-dir exp confirm
                confirm_up = True
            if dir_ == -1 and t_min > l_min and t_max > l_max:
                confirm_down = True
            # if have confirm change dir
            if confirm_up or confirm_down:
                b_flag = True
            if (confirm_up and t_max > l_max) or (confirm_down and t_min < l_min):
                c_flag = True
    
            if t_dir != dir_ and b_flag and not c_flag:
                exp_flag = True
            if exp_flag:
                macd_dir[l_ind:b_ind] = -macd_dir[l_ind:b_ind]
                #print(ind)
                break
    
            # update max min lmax, lmin
            max_dir[ind] = t_max
            min_dir[ind] = t_min
            max_last[ind] = l_max
            min_last[ind] = l_min
    
            # 保留上一部dir
            t_dir = dir_
    
            # all exp has refine
            if ind == macd_dir.shape[0] - 1:
                break_flag = False
    return macd_dir,max_dir,min_dir,max_last,min_last


def get_signal_p2(df,n0=60):
    for i,n in enumerate(range(n0,len(df)+1)):
        s=df['closeIndex'][:n].copy()
        macd_dir0 = macd_direction().integal_macd(s).copy()
        macd_dir0,max_dir0,min_dir0,max_last0,min_last0=get_signal(s,macd_dir0)
        if i ==0:
            macd_dir = macd_dir0
            max_dir = max_dir0
            min_dir = min_dir0
            max_last = max_last0
            min_last = min_last0
        else:
            macd_dir =np.hstack((macd_dir,macd_dir0[-1]))
            max_dir = np.hstack((max_dir,max_dir0[-1]))
            min_dir = np.hstack((min_dir,min_dir0[-1]))
            max_last = np.hstack((max_last,max_last0[-1]))
            min_last = np.hstack((min_last,min_last0[-1]))
        if np.mod(n,1000) == 0:
            print('S76part2信号计算：complete %d-%d' % (i,len(df)))
    # %%
    s_df = df[['closeIndex']].copy()
    s_df['dir'] = macd_dir
    s_df['max_dir'] = max_dir
    s_df['min_dir'] = min_dir
    s_df['max_last'] = max_last
    s_df['min_last'] = min_last
   # %% 策略exp  closeIndex > max_last表示方向确认，此时 exp 不起作用。没有涉及未来信息。
    #sign_long = s_df.apply(lambda x: 1 if x['dir'] == 1 and x['closeIndex'] > x['max_last'] else 0, axis=1)
    #signal1
    if method_sel==1:
        sign_long = s_df.apply(lambda x: 1 if x['dir'] == 1 and x['closeIndex'] > x['max_last'] \
                           else -1 if x['dir'] == -1 and x['closeIndex'] < x['max_last'] else 0, axis=1)
    elif method_sel==2:
    #signal2
        s_df['ma50']=s_df.closeIndex.rolling(50).mean()
        s_df['ma10']=s_df.closeIndex.rolling(10).mean()
        sign_long = s_df.apply(lambda x: 1 if x['dir'] == 1 and x['closeIndex'] > x['max_last'] and x['ma10']>x['ma50'] \
                           else -1 if x['dir'] == -1 and x['closeIndex'] < x['max_last'] and x['ma10']<=x['ma50'] else 0, axis=1)
        
        
    return sign_long

def get_t0_p2(ticker,dtype):
    sql_tmp = 'select * from %s where ticker = "%s" and dtype = "%s" order by tradeDate'
    df0 = pd.read_sql(sql_tmp % (tn,ticker,dtype),eg)
    if len(df0)>0:
        t0 = df0.tradeDate.astype(str).max()
    else:
        t0 = t00
    return t0,df0


def show_figure_p2(df0,dtype,ticker,disp):
    if disp:
        plt.figure(figsize=(10,5))    
        (df0['CHGPct'] + 1).cumprod().plot(label='index')
        (df0['CHGPct'] * df0['sig'].shift(1) + 1).cumprod().plot(label='stra1')
        plt.title('%s-%s' % (dtype,ticker))
        plt.show()
def update_csi_index():
    obj_t.use('S76P2-国内指数')
    for ticker in ["000001","000300","000905"]:
        dtype = 'index_A'
        #df = pd.read_csv('./data/hs_300.csv', index_col=0)
        df = MktIdxdGet(ticker,t00,'2099-01-01','tradeDate,closeIndex,CHGPct')
        df.tradeDate=df.tradeDate.astype(str)
        #df = df[df.tradeDate>=t00]
        df.reset_index(drop=True,inplace=True)    
        t0,df0 = get_t0_p2(ticker,dtype)
        if t0==t00:
            n0=60
        else:
            n0=df[df.tradeDate==t0].index.values[0]
        
        sign_long = get_signal_p2(df,n0)
        signal = df[['tradeDate','CHGPct']]
        sign_long.name='sig'
        signal = pd.concat([signal,sign_long],axis=1)
        signal = signal.copy()
        signal['ticker'] = ticker
        signal['dtype'] = dtype
        signal  = signal[tn_var]
        result = signal[signal.tradeDate>t0]
        if len(result)>0:
            result.to_sql(tn,eg,if_exists='append',index=False,chunksize=3000,dtype=dtypedict1)
        #save
        df0=pd.concat([df0,result])
        df0.set_index('tradeDate',inplace=True)
        
        show_figure_p2(df0,dtype,ticker,disp)


def update_ISO_index():
    obj_t.use('S76P2-国际指数')
    for ticker in ["SPX","N225","HSI"]:
        dtype = 'iso'
        #df = pd.read_csv('./data/hs_300.csv', index_col=0)
        df = MktIdxdGet(ticker,t00,'2099-01-01','tradeDate,closeIndex,CHGPct')
        df.tradeDate=df.tradeDate.astype(str)
        #df = df[df.tradeDate>=t00]
        df.reset_index(drop=True,inplace=True)    
        t0,df0 = get_t0_p2(ticker,dtype)
        if t0==t00:
            n0=60
        else:
            n0=df[df.tradeDate==t0].index.values[0]
        
        sign_long = get_signal_p2(df,n0)
        signal = df[['tradeDate','CHGPct']]
        sign_long.name='sig'
        signal = pd.concat([signal,sign_long],axis=1)
        signal = signal.copy()
        signal['ticker'] = ticker
        signal['dtype'] = dtype
        signal  = signal[tn_var]
        result = signal[signal.tradeDate>t0]
        if len(result)>0:
            result.to_sql(tn,eg,if_exists='append',index=False,chunksize=3000,dtype=dtypedict1)
        #save
        df0=pd.concat([df0,result])
        df0.set_index('tradeDate',inplace=True)
        
        show_figure_p2(df0,dtype,ticker,disp)
def update_forex():
    obj_t.use('S76P2-外汇')
    dtype = 'forex'
    ticker_pool = ['EURUSD', 'GBPUSD', 'USDJPY', 'AUDUSD', 'NZDUSD', 'USDCAD',
                   'USDCHF', 'EURGBP', 'EURCHF', 'EURAUD', 'EURNZD', 'EURCAD', 'AUDJPY', 'NZDJPY',
                   'EURJPY', 'CHFJPY', 'GBPJPY', 'CADJPY', 'USDCNH', 'USDINR', 'USDTRY', 'USDRUB',
                   'USDZAR']
    eg_polygon = create_db('polygon')
    sql_tmp = 'select tradeDate,closePrice as closeIndex from forex_day where ticker = "%s" order by tradeDate'
    for ticker in ticker_pool:
        df = pd.read_sql(sql_tmp % ticker,eg_polygon)
        df.tradeDate = df.tradeDate.astype(str)
        df['CHGPct'] = df.closeIndex.pct_change()
        df.reset_index(drop=True,inplace=True)    
        t0,df0 = get_t0_p2(ticker,dtype)
        if t0 > df.tradeDate.max():
            continue
        if t0==t00:
            n0=60
        else:
            n0=df[df.tradeDate==t0].index.values[0]
        
        sign_long = get_signal_p2(df,n0)
        signal = df[['tradeDate','CHGPct']]
        sign_long.name='sig'
        signal = pd.concat([signal,sign_long],axis=1)
        signal = signal.copy()
        signal['ticker'] = ticker
        signal['dtype'] = dtype
        signal  = signal[tn_var]
        result = signal[signal.tradeDate>t0]
        if len(result)>0:
            result.to_sql(tn,eg,if_exists='append',index=False,chunksize=3000,dtype=dtypedict1)
        #save
        df0=pd.concat([df0,result])
        df0.set_index('tradeDate',inplace=True)
        
        show_figure_p2(df0,dtype,ticker,disp)
def update_csi_future():
    obj_t.use('S76P2-国内期货')
    #国内期货
    #涨跌幅，（收盘价-昨结算价）/昨结算价
    dtype='csifuture'
    engine = create_db('yuqerdata')
    sql_tmp = 'select distinct(concat(contractObject,"-",exchangeCD)) from yq_MktMFutdGet'
    info = pd.read_sql(sql_tmp,engine)
    info = info[info.columns[0]].tolist()
    sql_1 ='''select  ticker, tradeDate,chgPct as CHGPct from yq_MktMFutdGet 
            where tradeDate>="%s" and contractObject="%s" and exchangeCD="%s" and mainCon=1 order by tradeDate'''
    for sub_info in info:
        ticker,exchangeCD = sub_info.split('-')
        df = pd.read_sql(sql_1 % (t00,ticker,exchangeCD),engine)
        if len(df)>440:
            df.tradeDate = df.tradeDate.astype(str)
            df['CHGPct'][df.ticker!=df.ticker.shift(1)]=0
            df['closeIndex'] = (1+df['CHGPct']).cumprod()
            df.reset_index(drop=True,inplace=True)    
            t0,df0 = get_t0_p2(ticker,dtype)
            if t0==t00:
                n0=60
            else:
                n0=df[df.tradeDate==t0].index.values[0]
            
            sign_long = get_signal_p2(df,n0)
            signal = df[['tradeDate','CHGPct']]
            sign_long.name='sig'
            signal = pd.concat([signal,sign_long],axis=1)
            signal = signal.copy()
            signal['ticker'] = ticker
            signal['dtype'] = dtype
            signal  = signal[tn_var]
            result = signal[signal.tradeDate>t0]
            if len(result)>0:
                result.to_sql(tn,eg,if_exists='append',index=False,chunksize=3000,dtype=dtypedict1)
            #save
            df0=pd.concat([df0,result])
            df0.set_index('tradeDate',inplace=True)            
            show_figure_p2(df0,dtype,ticker,disp)
#long-only
def bac_sim(df,ticker,dtype,t0):
    sub_fee1 = fee1[dtype]/10000
    sub_fee2 = fee2[dtype]/10000
    if len(df)>440:
        df.tradeDate = df.tradeDate.astype(str)
        df['CHGPct'][df.ticker!=df.ticker.shift(1)]=0
        df['closeIndex'] = (1+df['CHGPct']).cumprod()
        df.reset_index(drop=True,inplace=True)
        #t0,df0 = get_t0_p2(ticker,dtype)
        if t0==t00:
            n0=60
        else:
            n0=df[df.tradeDate==t0].index.values[0]
        
        sign_long = get_signal_p2(df,n0)
        signal = df[['tradeDate','CHGPct']]
        sign_long.name='sig'
        signal = pd.concat([signal,sign_long],axis=1)
        signal = signal.copy()
        signal['ticker'] = ticker
        signal['dtype'] = dtype
        signal  = signal[tn_var]
        signal.sig[signal.sig<0] = 0
        signal['sig1'] = signal['sig'].shift(1)
        signal['sig1'].fillna(0,inplace=True)
        signal['fee'] = 0
        signal.fee[(signal.sig1==1)&(signal.sig1.shift(1)==0)] = sub_fee1 #建仓
        signal.fee[(signal.sig1==0)&(signal.sig1.shift(1) == 1)] = sub_fee2 #清仓
        result = signal[signal.tradeDate>t0]    
        return result
    else:
        return pd.DataFrame()        

def bac_BX(data):
    p0,dtype,t0=data
    ticker,df = p0
    df.reset_index(drop = True,inplace=True)
    print('bac Test S76 %s-%s' % (dtype,ticker))
    return bac_sim(df,ticker,dtype,t0)

            
def update_csi_stock():
    obj_t.use('S76P2-国内股票')
    #必须实现并行，否者无法
    #国内期货
    dtype='csi_stock' 
    tn = 's76_%s' % dtype
    #engine = create_db('yuqerdata')    
    ticker_pool = get_symbol_A()
    ticker_pool = [i for i in ticker_pool if len(i) ==6]
    if table_in_database(dn,tn):
        sql_tmp = 'select ticker,max(tradeDate) as t from %s group by ticker'
        t_ini = pd.read_sql(sql_tmp % tn,eg)
        ## 需要补充
        if len(t_ini) == 0:
            t_ini = [t00] * len(ticker_pool)
        else:
            tmp = dict(zip(t_ini.ticker,t_ini.t.astype(str)))
            t_ini = []
            for i in ticker_pool:
                if i in tmp.keys():
                    t_ini.append(tmp[i])
                else:
                    t_ini.append(t00)
    else:
        ini_table_76p2(tn)
        t_ini = [t00] * len(ticker_pool)
    tt = get_table_date('yq_MktEqudAdjAfGet',engine,'tradeDate')
    if tt > max(t_ini):
        t_info = dict(zip(ticker_pool,t_ini))
        sub_str0 = 'ticker, tradeDate,closePrice as closeIndex,closePrice/preClosePrice-1 as CHGPct'
        if BX_sel:
            p1 = get_MktEqudAdjAfGet(t00,'2099-01-01',sub_str0)
            p1 = p1.groupby('ticker')
            sub_ticker = list(p1.groups.keys())
            p2 = [dtype] * len(sub_ticker)
            p3 = []
            for sub_sub_ticker in sub_ticker:
                if sub_sub_ticker in t_info.keys():
                    p3.append(t_info[sub_sub_ticker])
            #para
            pool = multiprocessing.Pool(num_core2)
            Y = pool.map(bac_BX, zip(p1,p2,p3))
            pool.close()
            pool.join()
        else:
            Y =[]
            for ticker in tqdm(ticker_pool):
                df = get_MktEqudAdjAfGet_update(ticker,t00,'2099-01-01',field=sub_str0)
                y = bac_sim(df,ticker,dtype,t_info[ticker])
                Y.append(y)
        #save_pickle('tmp.pkl',Y)
        Y=pd.concat(Y)
        if len(Y)>0:
            
            Y.to_sql(tn,eg,if_exists='append',index=False,chunksize=3000)


def get_single_s76(data):
    ticker,begin,end,dtype,t0=data
    sql_tmp = '''select ticker,tradeDate,closePrice from %s where ticker = "%s" and index_id="%s" 
        and tradeDate>="%s" and tradeDate<="%s" and closePrice is not null order by tradeDate'''
    df = pd.read_sql(sql_tmp % (tn_tdx,ticker,dtype,begin,end),eg_pro)
    df.drop_duplicates(subset=['ticker','tradeDate'],inplace=True)
    df.reset_index(drop=True,inplace=True)
    df.rename(columns={'closePrice':'closeIndex'},inplace=True)
    df['CHGPct'] = df.closeIndex.pct_change()
    df.CHGPct.fillna(0,inplace=True)
    y = bac_sim(df,ticker,dtype,t0)
    return y
    
    
            
def update_iso_com():
    time_start = datetime.datetime.now()
    for dtype in index_tdx:
        if dtype in dtype in ["HK",'hsce','hk_ggt']:
            continue
        tn = 's76_%s' % dtype.lower()
        begin = t00   
        end = datetime.datetime.strftime(datetime.datetime.now(),'%Y-%m-%d')
        if table_in_database('s37',tn):
            t0 = pd.read_sql('select tradeDate from %s where index0 = "%s" order by tradeDate desc limit 1' % (tn,dtype),engine37)
            if len(t0)>0:
                t0 = t0.tradeDate.astype(str)[0]
            else:
                t0 = begin
        else:
            #create_s43_tb_v2()
            t0 = begin
        tt = pd.read_sql('select tradeDate from main_index_s68 where index_id="%s" order by tradeDate desc limit 1' % dtype,eg_pro)
        tt = tt[tt.columns[0]].astype(str).values[0]
        if tt> t0:
            print('dealing with %s' % dtype)
            hs_300_pool=get_ticker_pool(dtype)
                
            T_symbols = len(hs_300_pool)
            if BX_sel:
                p1 = hs_300_pool
                p2 = T_symbols*[begin]
                p3 = T_symbols*[end]
                p4 = T_symbols*[dtype]
                p5 = T_symbols*[t0]
                pool = multiprocessing.Pool(num_core2)
                back_dfs = pool.map(get_single_s76, zip(p1,p2,p3,p4,p5))
                pool.close()
                pool.join() 
            else:    
                back_dfs = []
                for ticker in tqdm(hs_300_pool):
                    df_back = get_single_s76([ticker,begin,end,dtype,t0])
                    back_dfs.append(df_back)
            #sig0 = [i[1] for i in back_dfs]
            #back_dfs = [i[0] for i in back_dfs]
            #记录了信号的初始位置
            #sig0 = pd.concat(sig0)   
            #save_pickle('tmp.pkl',back_dfs)
            back_df = pd.concat(back_dfs)
            
            back_df=back_df[back_df.tradeDate.astype(str)>t0]
            #sig0=sig0[sig0.tradeDate.astype(str)>t0]
            if len(back_df)>0:
                back_df['index0'] = dtype
                #sig0['index0'] = dtype
                back_df.to_sql(tn,engine37,if_exists='append',index=False,chunksize=3000)
                #sig0.to_sql(tn_sig,engine37,if_exists='append',index=False,chunksize=3000)                
            time_end = datetime.datetime.now()
            print('Time used %s' % (time_end-time_start))
        else:
            print('S76 %s 数据已经是最新，无需更新' % dtype)
     
     
        
if __name__ == "__main__":
    # %%
    update_csi_index()
    update_ISO_index()
    update_forex()
    update_csi_future()

    update_csi_stock()
    update_iso_com()
    
    obj_t.use('S76 complete')
    
    
    