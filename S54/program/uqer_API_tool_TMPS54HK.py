# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 18:39:02 2020
S56历史数据临时程序
20201007 升级港股数据
@author: adair2019
"""
import uqer
from uqer import DataAPI
client = uqer.Client(token='8b1df403d1bfb39b588d35a43e8526f454d779351d0b439b9ad7b299cd61d9df')
import os
from yq_toolsS45 import get_file_name
from yq_toolsS45 import time_use_tool 
import datetime
obj_t= time_use_tool()
import pandas as pd
from multiprocessing.dummy import Pool as ThreadPool
import multiprocessing


#from yq_toolsSFZ import engine
server_sel = True
if server_sel:
    from yq_toolsSFZ import engine
    from yq_toolsSFZ import pn as datadir
else:
    from yq_toolsS45 import engine
    from yq_toolsS45 import pn as datadir

num_core = min(multiprocessing.cpu_count()*10,20)

#datadir='dataset_uqer'
if not os.path.exists(datadir):
    os.makedirs(datadir)


_,fn_exist=get_file_name(datadir,'.csv')

def get_ini_data(tn,var_name,db=engine):
    sql_str = 'select %s from %s order by %s desc limit 1' % (var_name,tn,var_name)
    t0 = pd.read_sql(sql_str,db)
    return t0[var_name].astype(str).values[0]
def list_files():
    return fn_exist

def save_data_adair(fn_d1,x,fn_d2=None):
    if fn_d2 is None:
        fn_d2 = fn_d1
    fn1_d1 =os.path.join(datadir, '%s.csv' % fn_d1)
    #fn2_d1 = '%s.zip' % fn_d2
    x.to_csv(fn1_d1,index=False)

def tref_split(tref,r=20):
    t0_1=[]
    tt_1=[]
    i=0    
    T = len(tref)
    while i <T-1:
        j=i+r-1
        if j>T-1:
            j=T-1
        t0_1.append(tref[i])
        tt_1.append(tref[j])
        i=i+r       
    return t0_1,tt_1
    
def get_symbol_adair():
    x=DataAPI.EquGet(secID=u"",ticker=u"",equTypeCD=u"A",listStatusCD=u"",field=u"",pandas="1")
    y=x['ticker']
    z=''
    z1=[]
    k=0;
    for i in y:
        if len(i)==6:
            if i[0]=='3' or i[0]=='0' or i[0:2]=='60':
                z1.append(i)
                if k==0:
                    z=i
                else:
                    z=z+','+i
                k = k +1
    return z1,z

fund_secID = DataAPI.FundGet(secID=u"",ticker=u"",etfLof=u"",listStatusCd=u"",
                        category=['E','H','B','SB','M','O'],idxID=u"",idxTicker=u"",
                        operationMode=u"",beginDate=u"",endDate=u"",status="",field='secID',pandas="1")
fund_secID=fund_secID.secID.unique().tolist()


z1,z=get_symbol_adair()

def get_tradingdate_adair(tt):
    x=DataAPI.TradeCalGet(exchangeCD=u"XSHG",beginDate=u"20000101",endDate=tt,field=u"calendarDate,isOpen",pandas="1")
    t=x.calendarDate[x.isOpen==1].values
    return t


tickerHK0=DataAPI.HKEquGet(secID=u"",ticker=u"",listStatusCD=u"",ListSectorCD=[1,2],equTypeCD=u"",connect=u"",field=u"",pandas="1")
tickerHK=tickerHK0.ticker.tolist()
trefHK0=DataAPI.TradeCalGet(exchangeCD=u"XHKG",beginDate=u"19900101",endDate=u"20300101",isOpen=u"",field=u"",pandas="1")
trefHK = trefHK0[trefHK0.isOpen==1].calendarDate.tolist()



class uq_methods:
    ####20
    def get_rontie_20d(self,t0,tt,para_sel=True):
        #1 S53 每日行情填充数据
        def get_MktEqudGet0S53():
            #股票日行情 半年
            key_str = '正股 filling 日数据'
            fn_d2 = 'MktEqudGet0S53_%s' % tt
            fn1_d2 = '%s.csv' % fn_d2
            fn2_d2= fn1_d2
            if fn2_d2 not in list_files():
                z1,z=get_symbol_adair()
                x=DataAPI.MktEqudGet(secID=u"",ticker=z1,tradeDate=u"",beginDate=t0,endDate=tt,isOpen=0,field=u"",pandas="1") 
                x.drop(['secID','secShortName','exchangeCD','vwap','isOpen'],axis=1, inplace=True)
                save_data_adair(fn_d2,x,'get_MktEqudGet0S53')
                print('%s已经更新到%s' % (key_str,fn2_d2)) 
            else:
                print('%s已经存在，%s已经更新，未执行' % (fn2_d2,key_str)) 
        #2
        def NewsHeatIndexNewGet_adair():
            if tt <='2015-04-03':
                return        
            info = '新闻热度指数（新版，不包括当天）'
            fn_d1 = 'NewsHeatIndexNewGetS55FZ_%s' % tt
            fn2_d2=  '%s.csv' % fn_d1
            if fn2_d2 in list_files():
                return
            x=DataAPI.NewsHeatIndexNewGet(beginDate=t0,endDate=tt,secID="",exchangeCD=u"",ticker=u"",secShortName=u"",field=u"",pandas="1")
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn_d1))  
        #3 新闻情感指数（新版，按天统计） 最好20-40天间隔
        def NewsSentiIndexGet_adair():
            if tt <='2015-04-03':
                return
            info = '新闻情感指数（新版，按天统计）'
            fn_d1 = 'NewsSentiIndexGetS55FZ_%s' % tt
            fn2_d2=  '%s.csv' % fn_d1
            if fn2_d2 in list_files():
                return
            
            x=DataAPI.NewsSentiIndexGet(beginDate=t0,endDate=tt,secID="",exchangeCD=u"",ticker=u"",
                                        secShortName=u"",field=u"",pandas="1")
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn_d1)) 
        #4 一致预期数据表(申万行业) 20-40
        def ResConInduSwGet_adair():
            if tt < '2014-01-02':
                return
            info = '一致预期数据表(申万行业)'
            fn_d1 = 'ResConInduSwGet18_%s' % tt
            fn2_d2=  '%s.csv' % fn_d1
            if fn2_d2 in list_files():
                return
            x= DataAPI.ResConInduSwGet(beginDate=t0,endDate=tt,secCode=u"",secName="",field="",pandas="1") 
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn_d1)) 
        #10-35
        #5 个股一致预期衍生数据表
        def ResConSecDerivativeGet_adair():
            if tt<'2010-01-04':
                return
            
            info = '个股一致预期衍生数据表'
            fn_d1 = 'ResConSecDerivativeGet18_%s' % tt
            fn2_d2=  '%s.csv' % fn_d1
            if fn2_d2 in list_files():
                return
            x=DataAPI.ResConSecDerivativeGet(secCode=z1,secName=u"",endDate=u"20200920",beginDate=u"20200910",field=u"",pandas="1")
            
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn_d1))  
            
        #max(41) 20比较合适
        #6 获取一致预期目标价与评级表
        def ResConTarpriScoreGet_adair():
            if tt<'20050104':
                return
            info='获取一致预期目标价与评级表'
            fn_d1='ResConTarpriScoreGet18_%s' % tt
            fn2_d2=  '%s.csv' % fn_d1
            if fn2_d2 in list_files():
                return
            
            x=DataAPI.ResConTarpriScoreGet(secCode=u"",endDate=tt,beginDate=t0,field="",pandas="1")   
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn_d1))
        #7        
        def get_MktEquFlowGet():
            if tt<'20090105':
                return
            info='个股日资金流向'
            fn_d1='MktEquFlowGetS56%s' % tt
            fn2_d2=  '%s.csv' % fn_d1
            if fn2_d2 in list_files():
                return
            
            x=DataAPI.MktEquFlowGet(tradeDate=u"",secID=u"",ticker=u"",
                                    beginDate=t0,endDate=tt,field=u"",pandas="1")  
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn_d1))
        #8
        def get_FstDetailGet():
            if tt<'20090102':
                return
            info='沪深融资融券每日交易明细信息'
            fn_d1='FstDetailGetS56%s' % tt
            fn2_d2=  '%s.csv' % fn_d1
            if fn2_d2 in list_files():
                return
            
            x=DataAPI.FstDetailGet(secID=u"",ticker=u"",beginDate=t0,endDate=tt,
                                   field=u"",pandas="1")
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn_d1))
        #9
        def get_HKshszHoldGet():
            if tt<'20170317':
                return
            info='沪深港通持股记录'
            fn_d1='HKshszHoldGetS56%s' % tt
            fn2_d2=  '%s.csv' % fn_d1
            if fn2_d2 in list_files():
                return
            
            x=DataAPI.HKshszHoldGet(secID=u"",ticker=u"",tradeCD=["1","2","3","4","5"],
                        ticketCode=u"",partyName=u"",beginDate=t0,endDate=tt,field=u"",pandas="1")
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn_d1))
        #10 x = DataAPI.MktHKEqudGet(secID=u"",ticker=tickerHK,tradeDate=u"",beginDate=t0,endDate=tt,field=u"",pandas="1")
        def get_MktHKEqudGet():
            if tt<'20010102':
                return
            info='港股日行情'
            fn_d1='MktHKEqudGetS54%s' % tt
            fn2_d2=  '%s.csv' % fn_d1
            if fn2_d2 in list_files():
                return
            
            x=DataAPI.MktHKEqudGet(secID=u"",ticker=tickerHK,tradeDate=u"",
                                           beginDate=t0,endDate=tt,field=u"",pandas="1")
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn_d1))
            
        
        def do_update(ind):
            if ind==0:
                get_MktEqudGet0S53()
            elif ind==1:
                NewsHeatIndexNewGet_adair()
            elif ind==2:
                NewsSentiIndexGet_adair()
            elif ind==3:
                ResConInduSwGet_adair()
            elif ind==4:
                ResConSecDerivativeGet_adair()
            elif ind==5:
                ResConTarpriScoreGet_adair()
            elif ind==6:
                get_MktEquFlowGet()
            elif ind==7:
                get_FstDetailGet()
            elif ind==8:
                get_HKshszHoldGet()
            elif ind==9:
                get_MktHKEqudGet()
                
        #w_n = 9        
        #ind_p=range(w_n+1) 
        ind_p=[9]
        if para_sel:
            #多进程
            pool = ThreadPool(processes=len(ind_p))
            pool.map(do_update, ind_p)
            pool.close()
            pool.join()
        else:
            for i in ind_p:
                do_update(i)
        
    ####5
    def get_rontie_3d(self,t0,tt,para_sel=True):        
        def ResConSecDataGet_adair():
            if tt < '2005-01-04':
                return
            info = '一致预期个股数据表'
            fn_d1 = 'ResConSecDataGet18_%s' % tt
            fn2_d2=  '%s.csv' % fn_d1
            if fn2_d2 in list_files():
                return
            x= DataAPI.ResConSecDataGet(secCode=u"",endDate=tt,beginDate=t0,field="",pandas="1") 
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn_d1)) 
        #0-5
        #获取一致预期个股营业收入表
        def ResConSecIncomeGet_adair():
            if tt<'20050104':
                return
            info='获取一致预期个股营业收入表'
            fn_d1='ResConSecIncomeGetS18_%s' % tt
            fn2_d2=  '%s.csv' % fn_d1
            if fn2_d2 in list_files():
                return
            x=DataAPI.ResConSecIncomeGet(secCode=u"",endDate=tt,beginDate=t0,field="",pandas="1")        
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn_d1)) 
        
        #基金历史净值(货币型,短期理财债券型除外)
        def get_FundNavGet():
            if tt< '20050104':
                return
            fn_d2= 'FundNavGet_S51_%s' % tt
            #fn1_d1 = '%s.csv' % fn_d1
            #fn2_d1 = '%s.zip' % fn_d1
            info = '基金历史净值(货币型,短期理财债券型除外)'
            if '%s.csv' % fn_d2 in list_files():
                return
            else:
                x=DataAPI.FundNavGet(secID=u"",ticker=u"",dataDate=u"",beginDate=t0,endDate=tt,
                                     partyID="",partyShortName="",field=u"",pandas="1")
                save_data_adair(fn_d2,x,'get_FundNavGet')
                print('%s已经更新到%s' % (info,fn_d2))  
                return x
        
        def get_MktIdxdGet():
            fn_d2= 'indicator_data%s' % tt
            info = '指数日线'
            if '%s.csv' % fn_d2 in list_files():
                return
            else: 
                x=DataAPI.MktIdxdGet(indexID=u"",ticker=u"",tradeDate=u"",beginDate=t0,endDate=tt,exchangeCD=u"",field=u"",pandas="1")
                save_data_adair(fn_d2,x,'indicator_data')
                print('%s已经更新到%s' % (info,fn_d2))  
                return x
        def do_update(ind):
            if ind==0:
                ResConSecDataGet_adair()
            elif ind==1:
                ResConSecIncomeGet_adair()
            elif ind==2:
                get_FundNavGet()
            elif ind==3:
                get_MktIdxdGet()
        w_n = 3        
        ind_p=range(w_n+1)   
        if para_sel:
            #多进程
            pool = ThreadPool(processes=w_n+1)
            pool.map(do_update, ind_p)
            pool.close()
            pool.join()
        else:
            for i in ind_p:
                do_update(i)
    
    ### 1天
    def get_rontie_1d(self,t0,tt,para_sel=True):
        def get_MktCmeFutdGet():
            if tt< '20160912':
                return
            fn_d2 = 'MktCmeFutdGet_S50%s' % tt
            if '%s.csv' % fn_d2 in list_files():
                return
            field = """ticker,tradeDate,deliYear,deliMonth,contractObject,preSettlePrice,
                preOpenInt,openPrice,highestPrice,highestPriceSide,lowestPrice,lowestPriceSide,closePrice,
                closePriceSide,settlePrice,chg,turnoverVol"""
                
            #sub_key = field.replace('\n','').strip().split(',')
            #sub_key = [i.strip() for i in sub_key]
            x =DataAPI.MktCmeFutdGet(ticker=u"",tradeDate=tt,beginDate=u"",endDate=u"",contractObject=u"",field=field,pandas="1")
            save_data_adair(fn_d2,x,'get_MktCmeFutdGet')
            print('CME期货日行情%s' % fn_d2)
            
        def get_FundAssetsGet():
            sub_tt = t0
            fn_d1= 'yq_FundAssetsGet_S51_%s' % sub_tt
            fn1_d1 = '%s.csv' % fn_d1
            fn2_d1 = fn1_d1
            info = '基金资产配置'
            if fn2_d1 not in list_files():
                x=DataAPI.FundAssetsGet(secID=u"",ticker=u"",reportDate=u"",updateTime=u"",
                                        beginDate=sub_tt,endDate=sub_tt,field=u"",pandas="1") 
                save_data_adair(fn_d1,x,'get_FundAssetsGet')
                print('%s已经更新到%s' % (info,fn2_d1))  
                return x
            else:
                print('%s已经存在，%s数据已经更新，未执行' % (fn2_d1,info)) 
                return None
        def get_MktStockFactorsOneDayProGet_full():
            if t0<'20070104':
                return
            sub_tt=t0
            #sub_t0 = datetime.datetime.strptime(t0_0, '%Y-%m-%d').strftime('%Y-%m-%d')
            #sub_tt = datetime.datetime.strptime(tt, '%Y%m%d').strftime('%Y-%m-%d')
            fn_d1= 'MktStockFactorsOneDayProGet%s' % sub_tt
            fn1_d1 = '%s.csv' % fn_d1
            fn2_d1 = fn1_d1
            info = 'MktStockFactorsOneDayProGet_full'
            key_str = u""
            if fn2_d1 not in list_files():
                sub_t=sub_tt
                sub_t = sub_t.replace('-','')
                x=DataAPI.MktStockFactorsOneDayProGet(tradeDate=sub_t,secID=u"",ticker=z1,field=key_str,pandas="1")
                save_data_adair(fn_d1,x,'MktStockFactorsOneDayProGet')
                print('%s已经更新到%s' % (info,fn2_d1))  
                return x
            else:
                print('%s已经存在，未执行' % (fn2_d1)) 
                return None
        def get_MktStockFactorsOneDayGet_full():
            sub_tt=t0
            #sub_t0 = datetime.datetime.strptime(t0_0, '%Y-%m-%d').strftime('%Y-%m-%d')
            #sub_tt = datetime.datetime.strptime(tt, '%Y%m%d').strftime('%Y-%m-%d')
            fn_d1= 'MktStockFactorsOneDayGet_re_%s' % sub_tt
            fn1_d1 = '%s.csv' % fn_d1
            fn2_d1 = fn1_d1
            info = 'MktStockFactorsOneDayGet_full'
            key_str = u""
            if fn2_d1 not in list_files():
                sub_t=sub_tt
                sub_t = sub_t.replace('-','')
                x=DataAPI.MktStockFactorsOneDayGet(tradeDate=sub_t,secID=u"",ticker=z1,field=key_str,pandas="1")
                save_data_adair(fn_d1,x,'get_MktStockFactorsOneDayGet_full')
                print('%s已经更新到%s' % (info,fn2_d1))  
                return x
            else:
                print('%s已经存在，未执行' % (fn2_d1)) 
                return None
        def get_FundHoldingsGet():
            fn_d2 = 'FundHoldingsGet_S51%s' % tt
            if '%s.csv' % fn_d2 in list_files():
                return
            x=DataAPI.FundHoldingsGet(secID=fund_secID,ticker=u"",reportDate=tt,beginDate=u"",endDate=u"",secType="",
                                      holdingTicker="",holdingSecID="",field=u"",pandas="1")
            save_data_adair(fn_d2,x)
            print('基金持仓明细%s' % fn_d2)
        def do_update(ind):    
            if ind==0:
                get_MktCmeFutdGet()
            elif ind==1:
                get_FundAssetsGet()
            elif ind==2:
                get_MktStockFactorsOneDayProGet_full()
            elif ind==3:
                get_MktStockFactorsOneDayGet_full()
            elif ind==4:
                #15基金持仓明细
                get_FundHoldingsGet()
        w_n = 4        
        ind_p=range(w_n+1)   
        if para_sel:
            #多进程
            pool = ThreadPool(processes=w_n+1)
            pool.map(do_update, ind_p)
            pool.close()
            pool.join()
        else:
            for i in ind_p:
                do_update(i)
    
    ###全部 1次
    #1
    def get_rontie_0(self,t0,tt,para_sel=True):
        def EcoDataProGet_update():
            fn_d2 = 'EcoDataProGet_S53_%s' % t0
            fn2_d2=  '%s.csv' % fn_d2
            if fn2_d2 in list_files():
                return
            x=DataAPI.EcoDataProGet(indicID=['1020000004', '1030000011', '1040000050', '1040000702',
                   '1070000007', '1070000009', '1090001390', '1090001558'], beginDate=t0, endDate=tt, field=u"")
            save_data_adair(fn_d2,x)
            print('S53宏观数据已经更新 %s' % tt) 
            
        #2股票基本信息
        def get_symbol_basic_info():     
            fn_d2 = 'EquGet%s' % tt
            fn2_d2=  '%s.csv' % fn_d2
            if fn2_d2 not in list_files():
                x = DataAPI.EquGet(secID=u"",ticker=u"",equTypeCD=u"A",listStatusCD=u"",
                           exchangeCD="",ListSectorCD=u"",field=u"",pandas="1")   
                save_data_adair(fn_d2,x,'EquGet')
                print('股票基本信息已经更新到%s' % fn_d2)
            else:
                print('%s已经存在，股票基本信息已经更新，未执行' % fn2_d2) 
                
        #3 S53 上市公司特殊状态
        def EquInstSstateGet_adair():
            fn_d2 = 'EquInstSstateGet%s' % tt
            fn2_d2=  '%s.csv' % fn_d2
            if fn2_d2 in list_files():
                return
            x=DataAPI.EquInstSstateGet(secID=u"",ticker=z1,beginDate=t0,endDate=tt,field=u"",pandas="1")
            save_data_adair(fn_d2,x)
            print('已经更新到%s' % (fn_d2))
            
        #4 S51 基金基本信息
        def get_FundGet():
            fn_d2 = 'FundGet_S51%s' % tt
            fn2_d2=  '%s.csv' % fn_d2
            if fn2_d2 in list_files():
                return
            
            fields = ['secID', 'ticker', 'secShortName', 'tradeAbbrName', 'category',
                   'operationMode', 'indexFund', 'etfLof', 'isQdii', 'isFof', 'isGuarFund',
                   'guarPeriod', 'guarRatio', 'exchangeCd', 'listStatusCd', 'managerName',
                   'status', 'establishDate', 'listDate', 'delistDate', 'expireDate',
                   'managementCompany', 'managementFullName', 'custodian',
                   'custodianFullName',  'perfBenchmark',
                   'circulationShares', 'isClass', 'idxID', 'idxTicker', 'idxShortName',
                   'managementShortName', 'custodianShortName']
            t=['E','H','B','SB','M','O']
            y = DataAPI.FundGet(secID=u"",ticker=u"",etfLof=u"",listStatusCd=u"",
                                category=t,idxID=u"",idxTicker=u"",operationMode=u"",
                                beginDate=u"",endDate=u"",status="",field=fields,pandas="1")
            save_data_adair(fn_d2,y)
            print('已经更新到%s' % (fn_d2))
            
        #5 ETF基金申赎清单基本信息
        def FundETFPRListGet_adair():
            fn_d2 = 'FundETFPRListGet_adair%s' % tt
            fn2_d2=  '%s.csv' % fn_d2
            if fn2_d2 in list_files():
                return
            x = pd.DataFrame()
            for sub_code in ['510050','510300','510500']:
                sub_x=DataAPI.FundETFPRListGet(secID=u"",ticker=sub_code,
                                               beginDate=t0,endDate=tt,field=u"",pandas="1")
                x = x.append(sub_x)
            save_data_adair(fn_d2,x)
            print('已经更新到%s' % (fn_d2))
            
        #6 期货合约信息 
        def FutuGet_adair():
            info = '期货合约信息'
            fn_d2 = 'FutuGet_adair%s' % tt
            fn2_d2=  '%s.csv' % fn_d2
            if fn2_d2 in list_files():
                return
            x=DataAPI.FutuGet(secID=u"",ticker=u"",exchangeCD=u"",contractStatus="",contractObject=u"",field=u"",pandas="1")
            x.drop(['deliGrade','deliPriceMethod','settPriceMethod'],axis=1,inplace=True)
            save_data_adair(fn_d2,x)
            print('%s已经更新到%s' % (info,fn_d2))
            
        #7 申万获取行业
        def get_industry_data_adair():    
            fn_d2 = 'EquIndustryGet%s' % tt
            fn2_d2=  '%s.csv' % fn_d2
            if fn2_d2 not in list_files():
                x =DataAPI.EquIndustryGet(secID=u"",ticker=u"",
                                      industryVersionCD=u"010303",industry=u"",industryID=u"",industryID1=u"",
                                      industryID2=u"",
                                      industryID3=u"",intoDate=u"",equTypeID=u"",field=u"",pandas="1")  
                save_data_adair(fn_d2,x)
                print('申万获取行业信息已经更新到%s' % fn_d2)
            else:
                print('%s已经存在，申万获取行业信息已经更新，未执行' % fn2_d2) 
        #8 S49 申万行业回填（含科创板）    
        def get_MdSwBackGet():
            fn_d2 = 'MdSwBackGet_data%s' % tt
            fn2_d2=  '%s.csv' % fn_d2
            if fn2_d2 in list_files():
                return
            x=DataAPI.MdSwBackGet(secID=u"",ticker=u"",intoDate=u"",outDate=u"",field=u"",pandas="1")
            save_data_adair(fn_d2,x)
            print('申万行业回填（含科创板）%s' % fn_d2)
        #9  交易日数据已更新  
        def get_TradeCalGet():
            fn_d2 = 'yuqer_cal%s' % tt
            fn2_d2=  '%s.csv' % fn_d2
            if fn2_d2 in list_files():
                return
            x=DataAPI.TradeCalGet(exchangeCD=u"XSHG,XSHE",beginDate=u"",endDate=tt,isOpen=u"",field=u"",pandas="1")
            save_data_adair(fn_d2,x)
            print('交易日数据已更新%s' % fn_d2)
    
        def do_update(ind):
            if ind==0:
                EcoDataProGet_update()
            elif ind==1:
                get_symbol_basic_info()
            elif ind==2:
                EquInstSstateGet_adair()
            elif ind==3:
                get_FundGet()
            elif ind==4:
                FundETFPRListGet_adair()
            elif ind==5:
                FutuGet_adair()
            elif ind==6:
                get_industry_data_adair()
            elif ind==7:
                get_MdSwBackGet()
        w_n = 9        
        ind_p=range(w_n+1)   
        if para_sel:
            #多进程
            pool = ThreadPool(processes=w_n+1)
            pool.map(do_update, ind_p)
            pool.close()
            pool.join()
        else:
            for i in ind_p:
                do_update(i)

        
    ######可以年度更新的数据
    def get_rontie_year(self,t0,tt,para_sel=True):
        #指数估值
        def MktIdxdEvalGet_adair():
            info='指数估值'
            fn_d2 = 'MktIdxdEvalGet_adair%s' % tt
            if '%s.csv' % fn_d2 in list_files():
                return
            x= DataAPI.MktIdxdEvalGet(secID=u"",ticker=u"000922",beginDate=t0,endDate=tt,PEType=u"",field=u"",pandas="1")
            save_data_adair(fn_d2,x,'MktIdxdEvalGet_adair')
            print('%s已经更新到%s' % (info,fn_d2))
            
        def IdxCloseWeightGet_adair(): 
            if tt<'20050101':
                return
            fn_d2 = 'IdxCloseWeightGet%s' % tt
            fn2_d2=  '%s.csv' % fn_d2
            if fn2_d2 in list_files():
                return
            tickercode = ['000001','000002','000003','000004','000005','000006','000007','000008','000009','000010','000011','000012','000013','000015',
                          '000016','000020','000090','000132','000133','000300','000852','000902','000903','000904','000905','000906','000907','000922',
                          '399001','399002','399004','399005','399006','399007','399008','399009','399010','399011','399012','399013','399015','399107',
                          '399108','399301','399302','399306','399307','399324','399330','399333','399400','399401','399649','000985']
            x=[]
            for i,sub_ticker in enumerate(tickercode):
                sub_x=DataAPI.IdxCloseWeightGet(secID=u"",ticker=sub_ticker,beginDate=t0,endDate=tt,field=u"effDate,ticker,consTickerSymbol,weight",pandas="1")
                x.append(sub_x)
                print('%d-%s' % (i,sub_ticker))
            x=pd.concat(x)
            save_data_adair(fn_d2,x,'IdxCloseWeightGet')
            print('%s已经更新' % fn_d2) 
        #ST标记
        def update_st_data():   
            fn_d2 = 'st_data%s' % tt
            fn2_d2=  '%s.csv' % fn_d2
            if fn2_d2 not in list_files():
                x=DataAPI.SecSTGet(beginDate=t0,endDate=tt,secID=u"",ticker=z1,
                                   field=['ticker','tradeDate','STflg'],pandas="1")
                save_data_adair(fn_d2,x,'SecSTGet')
                print('st交易日数据已经更新到%s' % fn_d2)
            else:
                print('%s已经存在，交易日数据已经更新，未执行' % fn2_d2)  
        #按照年更新没有问题
        #2000-01-01
        def MktIborGet_adair():
            t0_0=t0
            sub_t0=t0_0
            sub_tt=tt
            fn_d1= 'MktIborGet_adair_%s' % sub_tt
            fn2_d2 = '%s.csv' % fn_d1
            if fn2_d2 in list_files():
                return
            info = '银行间同业拆借利率'
            
            ticker=['Hibor10M', 'Hibor11M', 'Hibor1D', 'Hibor1M', 'Hibor1W', 'Hibor1Y', 
                    'Hibor2M', 'Hibor2W', 'Hibor3M', 'Hibor4M', 'Hibor5M', 'Hibor6M', 'Hibor7M', 'Hibor8M', 'Hibor9M',
                 'Libor10M', 'Libor11M', 'Libor1D', 'Libor1M', 'Libor1W', 'Libor1Y', 
                 'Libor2M', 'Libor2W', 'Libor3M', 'Libor4M', 'Libor5M', 'Libor6M', 'Libor7M', 'Libor8M', 'Libor9M', 'Shibor1D',
                 'Shibor1D10D', 'Shibor1D20D', 'Shibor1D5D', 'Shibor1M', 'Shibor1M10D',
                 'Shibor1M20D', 'Shibor1M5D', 'Shibor1W', 'Shibor1W10D', 'Shibor1W20D', 'Shibor1W5D', 'Shibor1Y', 'Shibor1Y10D',
                 'Shibor1Y20D', 'Shibor1Y5D', 'Shibor2W', 'Shibor2W10D', 'Shibor2W20D', 
                 'Shibor2W5D', 'Shibor3M', 'Shibor3M10D', 'Shibor3M20D', 'Shibor3M5D', 'Shibor6M', 'Shibor6M10D', 'Shibor6M20D',
                 'Shibor6M5D', 'Shibor9M', 'Shibor9M10D', 'Shibor9M20D', 'Shibor9M5D', 
                 'Tibor10M', 'Tibor11M', 'Tibor12M', 'Tibor1M', 'Tibor1W', 'Tibor2M', 'Tibor3M', 'Tibor4M', 'Tibor5M', 'Tibor6M',
                 'Tibor7M', 'Tibor8M', 'Tibor9M']
            x=DataAPI.MktIborGet(secID=u"",ticker=ticker,tradeDate=u"",beginDate=sub_t0,endDate=sub_tt,currency=u"",field=u"",pandas="1")
            save_data_adair(fn_d1,x,'MktIborGet_adair')
            print('%s已经更新到%s' % (info,fn_d1))
        def do_update(ind):    
            if ind == 0:
                MktIdxdEvalGet_adair()
            elif ind==1:
                IdxCloseWeightGet_adair()
            elif ind==2:
                update_st_data()
            elif ind==3:
                MktIborGet_adair()
        w_n = 3        
        ind_p=range(w_n+1)   
        if para_sel:
            #多进程
            pool = ThreadPool(processes=w_n+1)
            pool.map(do_update, ind_p)
            pool.close()
            pool.join()
        else:
            for i in ind_p:
                do_update(i)
    
    ###半年
    def get_rontie_halfyear(self,t0,tt,para_sel=True):
        def MktMFutdGet_adair():
            if tt<'20010102':
                return
            info='期货主力、连续合约日行情'
            fn_d2 = 'MktMFutdGet_adair%s' % tt
            fn2_d2=  '%s.csv' % fn_d2
            if fn2_d2 in list_files():
                return
            x=DataAPI.MktMFutdGet(mainCon=u"",contractMark=u"",contractObject=u"",
                                  tradeDate=u"",startDate=t0,endDate=tt,field=u"",pandas="1")
            save_data_adair(fn_d2,x,'MktMFutdGet_adair')
            print('%s已经更新到%s' % (info,fn_d2))    
        #半年
        def FundETFConsGet_adair():
            if tt<'20160104':
                return
            fn_d2 = 'FundETFConsGet_adair %s' % tt
            fn2_d2=  '%s.csv' % fn_d2
            if fn2_d2 in list_files():
                return
            x = []
            for sub_code in ['510050','510300','510500']:
                sub_x=DataAPI.FundETFConsGet(secID=u"",ticker=sub_code,beginDate=t0,endDate=tt,field=u"",pandas="1")
                x.append(sub_x)
            x=pd.concat(x)
            save_data_adair(fn_d2,x,'FundETFConsGet_adair')
            print('已经更新到%s' % (fn_d2)) 
        #后复权月度行情
        def get_month_data():    
            fn_d2 = 'MktEqumAdjAfGet%s' % tt
            fn2_d2=  '%s.csv' % fn_d2
            if fn2_d2 not in list_files():
                x=DataAPI.MktEqumAdjAfGet(secID=u"",ticker=z1,monthEndDate=u"",beginDate=t0,
                                          endDate=tt,isOpen=u"",field=u"",pandas="1")
                save_data_adair(fn_d2,x,'MktEqumAdjAfGet')
                print('后复权月度行情已经更新到%s' % fn_d2)
            else:
                print('%s已经存在，交易日数据已经更新，未执行' % fn2_d2)   
        #指数月度行情
        def MktIdxmGet_adair():
            fn_d2 = 'MktIdxmGet%s' % tt
            if '%s.csv' % fn_d2 in list_files():
                return
            x=DataAPI.MktIdxmGet(beginDate=t0,endDate=tt,indexID=u"",ticker=u"",field=u"",pandas="1")
            save_data_adair(fn_d2,x,'MktIdxmGet')
        #后复权周度行情
        def get_ticker_week_data():    
            fn_d2 = 'MktEquwAdjAfGet%s' % tt
            fn2_d2=  '%s.csv' % fn_d2
            if fn2_d2 not in list_files():
                x=DataAPI.MktEquwAdjAfGet(secID=u"",ticker=z1,weekEndDate=u"",beginDate=t0,
                                          endDate=tt,isOpen=u"",field=u"",pandas="1")
                save_data_adair(fn_d2,x,'get_ticker_week_data')
                print('后复权周度行情已经更新到%s' % fn_d2)
            else:
                print('%s已经存在，交易日数据已经更新，未执行' % fn2_d2)
            #指数周度行情
        def get_index_week_data():    
            fn_d2 = 'MktIdxwGet_adair%s' % tt
            fn2_d2=  '%s.csv' % fn_d2
            if fn2_d2 not in list_files():
                x=DataAPI.MktIdxwGet(beginDate=t0,endDate=tt,indexID=u"",ticker=u"",field=u"",pandas="1")
                save_data_adair(fn_d2,x,'get_index_week_data')
                print('后复权周度行情已经更新到%s' % fn_d2)
            else:
                print('%s已经存在，交易日数据已经更新，未执行' % fn2_d2) 
                
        def do_update(ind):
            if ind==0:
                MktMFutdGet_adair()
            elif ind==1:
                FundETFConsGet_adair()
            elif ind==2:
                get_month_data()
            elif ind==3:
                MktIdxmGet_adair()
            elif ind==4:
                get_ticker_week_data()
            elif ind==5:
                get_index_week_data()
        w_n = 5        
        ind_p=range(w_n+1)   
        if para_sel:
            #多进程
            pool = ThreadPool(processes=w_n+1)
            pool.map(do_update, ind_p)
            pool.close()
            pool.join()
        else:
            for i in ind_p:
                do_update(i)
            
    ###10天肯定没有问题
    def get_rontie_10d(self,t0,tt,para_sel=True):
        def get_ticker_data():
            #x.loc[x.tradeDate>='2019-12-01']
            #2 股票日行情
            #股票日行情 半年
            fn_d2 = 'tickerday_data%s' % tt
            fn1_d2 = '%s.csv' % fn_d2
            fn2_d2= fn1_d2
            if fn2_d2 not in list_files():
                z1,z=get_symbol_adair()
                x=DataAPI.MktEqudGet(secID=u"",ticker=z1,tradeDate=u"",beginDate=t0,endDate=tt,isOpen=1,field=u"",pandas="1") 
                x.drop(['secID','secShortName','exchangeCD','vwap','isOpen'],axis=1, inplace=True)
    
                save_data_adair(fn_d2,x,'MktEqudGet')
                print('正股日数据已经更新到%s' % fn2_d2) 
            else:
                print('%s已经存在，正股日数据已经更新，未执行' % fn2_d2) 
      
        #沪深股票后复权行情
        def MktEqudAdjAfGet():
            fn_d2 = 'MktEqudAdjAfGet_data%s' % tt
            fn1_d2 = '%s.csv' % fn_d2
            fn2_d2= fn1_d2
            if fn2_d2 not in list_files():
                z1,z=get_symbol_adair()
                sub_x=DataAPI.MktEqudAdjAfGet(secID=u"",ticker=z1,tradeDate=u"",beginDate=t0,endDate=tt,
                                              isOpen=1,field=u"",pandas="1")
                save_data_adair(fn_d2,sub_x,'MktEqudAdjAfGet')
                print('后复权因子数据已经更新到%s' % fn2_d2) 
            else:
                #return None
                print('%s已经存在，后复权因子数据已经更新，未执行' % fn2_d2) 
        #S13
        #基金日行情
        def MktFunddget_adair():   
            info='基金日行情'
            fn_d2 = 'MktFunddget_adair%s' % tt
            if '%s.csv' % fn_d2 in list_files():
                return
            x=DataAPI.MktFunddGet(secID=u"",ticker=u"",tradeDate=u"",beginDate=t0,endDate=tt,field=u"",pandas="1")
            save_data_adair(fn_d2,x,'MktFunddget_adair')
            print('%s已经更新到%s' % (info,fn_d2))   
        #期货期货会员成交量排名 x = DataAPI.MktFutMTRGet(beginDate=sub_tref[0],endDate=sub_tref[-1],secID=u"",ticker=u"",field=u"",pandas="1")
        def MktFutMTRGet_adair():
            info='期货期货会员成交量排名'
            fn_d2 = 'MktFutMTRGet_adair%s' % tt
            if '%s.csv' % fn_d2 in list_files():
                return
            x= DataAPI.MktFutMTRGet(beginDate=t0,endDate=tt,secID=u"",ticker=u"",field=u"",pandas="1")
            save_data_adair(fn_d2,x,'MktFutMTRGet_adair')
            print('%s已经更新到%s' % (info,fn_d2))  
    
        #期货期货会员空头持仓排名 x = DataAPI.MktFutMSRGet(beginDate=sub_tref[0],endDate=sub_tref[-1],secID=u"",ticker=u"",field=u"",pandas="1")
        def MktFutMSRGet_adair():
            info='期货期货会员空头持仓排名'
            fn_d2 = 'MktFutMSRGet_adair%s' % tt
            if '%s.csv' % fn_d2 in list_files():
                return
            x= DataAPI.MktFutMSRGet(beginDate=t0,endDate=tt,secID=u"",ticker=u"",field=u"",pandas="1")
            save_data_adair(fn_d2,x,'MktFutMSRGet_adair')
            print('%s已经更新到%s' % (info,fn_d2))  
    
        #期货期货会员多头持仓排名 x = DataAPI.MktFutMLRGet(beginDate=sub_tref[0],endDate=sub_tref[-1],secID=u"",ticker=u"",field=u"",pandas="1")
        def MktFutMLRGet_adair():
            info='期货期货会员空头持仓排名'
            fn_d2 = 'MktFutMLRGet_adair%s' % tt
            if '%s.csv' % fn_d2 in list_files():
                return
            x= DataAPI.MktFutMLRGet(beginDate=t0,endDate=tt,secID=u"",ticker=u"",field=u"",pandas="1")
            save_data_adair(fn_d2,x,'MktFutMLRGet_adair')
            print('%s已经更新到%s' % (info,fn_d2))
        
        #期货仓单日报
        def MktFutWRdGet_adair():
            info='期货仓单日报'
            fn_d2 = 'MktFutWRdGet_adair%s' % tt
            if '%s.csv' % fn_d2 in list_files():
                return
            x= DataAPI.MktFutWRdGet(beginDate=t0,endDate=tt,contractObject=u"",exchangeCD=u"",field=u"",pandas="1")
            save_data_adair(fn_d2,x,'MktFutWRdGet_adair')
            print('%s已经更新到%s' % (info,fn_d2))    
        #可转债市场表现
        def MktConsBondPerfGet_adair():    
            fn_d1= 'MktConsBondPerfGet_adair%s' % tt
            fn2_d1 = '%s.csv' % fn_d1
            if fn2_d1 not in list_files():
                z1,z=get_symbol_adair()
                x=DataAPI.MktConsBondPerfGet(beginDate=t0,endDate=tt,secID=u"",tickerBond=u"",tickerEqu=z1,field=u"",pandas="1")
                x.drop(["reviseItem","triggerItem","triggerCondItem"],axis=1,inplace=True)
                save_data_adair(fn_d1,x,'MktConsBondPerfGet_adair')
                print('可转债日数据已经更新到%s' % fn2_d1)    
            else:
                print('%s已经存在，可转债日数据已经更新，未执行' % fn2_d1) 
        #S46 贴吧数据 
        def get_SocialDataGubaGet():
            fn_d2 = 'SocialDataGubaGet_data%s' % tt
            if '%s.csv' % fn_d2 in list_files():
                return
            x=DataAPI.SocialDataGubaGet(ticker=u"",beginDate=t0,endDate=tt,field=u"",pandas="1")
            save_data_adair(fn_d2,x,'get_SocialDataGubaGet')
            print('贴吧数据已更新%s' % fn_d2)
           
        #获取后复权填充数据
        def MktEqudAdjAfGetF0S53():    
            key_str0 = 'MktEqudAdjAfGetF0S53'
            key_str1 ='后复权因子填充数据'
            fn_d2 = '%s_%s' % (key_str0,tt)
            fn1_d2 = '%s.csv' % fn_d2
            fn2_d2= fn1_d2
            if fn2_d2 not in list_files():
                z1,z=get_symbol_adair()
                #t= get_tradingdate_adair(tt)
                x=DataAPI.MktEqudAdjAfGet(secID=u"",ticker=z1,tradeDate=u"",beginDate=t0,endDate=tt,
                                              isOpen=0,field=u"",pandas="1") 
                save_data_adair(fn_d2,x,'MktEqudAdjAfGetF0S53')
                #return x
                print('%s已经更新到%s' % (key_str1,fn2_d2)) 
            else:
                #return None
                print('%s已经存在，%s已经更新，未执行' % (fn2_d2,key_str1))
        def do_update(ind):
            #1
            if ind==0:
                get_ticker_data()
            elif ind==1:
                #2
                MktEqudAdjAfGet()
            elif ind==2:
                #3
                MktFunddget_adair()
            elif ind==3:
                #4
                MktFutMTRGet_adair()
            elif ind==4:
                #5
                MktFutMSRGet_adair()
            elif ind==5:
                #6
                MktFutMLRGet_adair()
            elif ind==6:
                #7
                MktFutWRdGet_adair()
            elif ind==7:
                #8
                MktConsBondPerfGet_adair()
            elif ind==8:
                #9
                get_SocialDataGubaGet()
            elif ind==9:
                #10
                MktEqudAdjAfGetF0S53()
        w_n = 9        
        ind_p=range(w_n+1)   
        if para_sel:
            #多进程
            pool = ThreadPool(processes=w_n+1)
            pool.map(do_update, ind_p)
            pool.close()
            pool.join()
        else:
            for i in ind_p:
                do_update(i) 
        
        
    #升级财务
    #组内并行
    def get_finance_data(self,t0,tt,para_sel=True):
        def get_FdmtISGet():
            info = '合并利润表'
            fn_d1= 'FdmtISGet_%s' % tt
            fn2_d2 = '%s.csv' % fn_d1            
            if fn2_d2 in list_files():
                return                    
            x=DataAPI.FdmtISGet(ticker=z1,secID=u"",reportType=u"",endDate=tt,beginDate=t0,publishDateEnd=u"",publishDateBegin=u"",
                      endDateRep="",beginDateRep="",beginYear="",endYear="",fiscalPeriod="",field=u"",pandas="1")
            
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn_d1))  
            
        #2 资产重组
        #x=DataAPI.EquRestructuringGet(secID=u"",ticker=u"000040",beginDate=u"20141231",endDate=u"",field=u"",pandas="1")
        def get_EquRestructuringGet():
            fn_d1= 'EquRestructuringGet_%s' % tt
            fn2_d2 = '%s.csv' % fn_d1            
            if fn2_d2  in list_files():
                return
            info = '资产重组表'
            
            f_str = [u'secID', u'ticker', u'secShortName', u'exchangeCD', u'publishDate',
               u'iniPublishDate', u'finPublishDate', u'program', u'isSucceed',
               u'restructuringType', u'underlyingType', u'underlyingVal',
               u'expenseVal', u'isRelevance', u'isMajorRes', u'payType',
               u'institNameB', u'relationShipB', u'institNameS', u'relationShipS',
               u'institNameSub', u'relationShipSub', u'institNameDeb',
               u'relationShipDeb', u'institNameCred', u'relationShipCred']
            x=DataAPI.EquRestructuringGet(secID=u"",ticker=z1,beginDate=t0,endDate=tt,field=f_str,pandas="1")
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn_d1))  
            
            
        #3 合并资产负债表 (Point in time) 
        #AR	float	应收账款
        #TCA	float	流动资产合计
        #TAssets	float	资产总计
        #NotesReceiv	float	应收票据
        #othReceiv	float	其他应收款
        #othCA	float	其他流动资产
        #TCL	float	流动负债合计
        #DataAPI.FdmtBSGet(ticker=u"688001",secID=u"",reportType=u"",endDate=u"",beginDate=u"",publishDateEnd=u"",
        #publishDateBegin=u"",endDateRep="",beginDateRep="",beginYear="",endYear="",fiscalPeriod="",field=u"",pandas="1")
        def get_FdmtBSGet():
            fn_d1= 'FdmtBSGet_%s' % tt
            fn2_d2 = '%s.csv' % fn_d1            
            if fn2_d2 in list_files():
                return
            fn2_d1 = fn_d1
            info = '合并资产负债表'
            x=DataAPI.FdmtBSGet(secID=u"",ticker=z1,beginDate=t0,endDate=tt,publishDateBegin=u"",field=u"",pandas="1")
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn2_d1))  
            
        # 4主营业务构成（基础数据）
        #grossMargin	float	毛利率
        #DataAPI.FdmtMainOperNGet(partyID="",secID=u"",ticker=u"688001",beginDate=u"20181231",endDate=u"",field=u"",pandas="1")
        def get_FdmtMainOperNGet():
            fn_d1= 'FdmtMainOperNGet_%s' % tt
            fn2_d2 = '%s.csv' % fn_d1            
            if fn2_d2 in list_files():
                return
            fn2_d1 = fn_d1
            info = '主营业务构成'
            x=DataAPI.FdmtMainOperNGet(secID=u"",ticker=z1,beginDate=t0,endDate=tt,field=u"",pandas="1")
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn2_d1))  
            
        #业绩快报 数据太少，一直没有用
        #DataAPI.FdmtEeGet(ticker=u"600000",secID=u"",reportType=u"",endDate=u"",beginDate=u"",publishDateEnd=u"",publishDateBegin=u"",field=u"",pandas="1")
        def get_FdmtEeGet_S26():
            fn_d1= 'FdmtEeGet_S26_%s' % tt
            fn2_d2 = '%s.csv' % fn_d1            
            if fn2_d2  in list_files():
                return
            fn2_d1 = fn_d1
            info = '业绩快报'
            x=DataAPI.FdmtEeGet(secID=u"",ticker=z1,beginDate=t0,endDate=tt,field=u"",pandas="1")
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn2_d1))  
            
        #DataAPI.FdmtDerPitGet(secID="",ticker=u"688002",beginDate="",endDate="",beginYear=u"",endYear=u"",reportType=u"",publishDateEnd=u"",publishDateBegin=u"",field=u"",pandas="1")
        #5 财务衍生数据 (Point in time)
        #nrProfitLoss 非经常性损益 , 直接取公告披露值
        def get_FdmtDerPitGet():
            fn_d1= 'FdmtDerPitGet_%s' % tt
            fn2_d2 = '%s.csv' % fn_d1            
            if fn2_d2 in list_files():
                return
            fn2_d1 = fn_d1
            info = '财务衍生数据 (Point in time)'
            #x=DataAPI.FdmtEeGet(secID=u"",ticker=z1,beginDate=sub_t0,endDate=sub_tt,field=u"",pandas="1")
            x=DataAPI.FdmtDerPitGet(secID="",ticker=z1,beginDate=t0,endDate=tt,beginYear=u"",endYear=u"",
                                    reportType=u"",publishDateEnd=u"",publishDateBegin=u"",field=u"",pandas="1")
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn2_d1))  
            
        #6 财务指标—营运能力 (Point in time)
        #[通联数据] - DataAPI.FdmtIndiTrnovrPitGet
        def get_FdmtIndiTrnovrPitGet():
            fn_d1= 'FdmtIndiTrnovrPitGet_%s' % tt
            fn2_d2 = '%s.csv' % fn_d1            
            if fn2_d2 in list_files():
                return
            fn2_d1 = fn_d1
            info = '财务指标-运营能力'
            x=DataAPI.FdmtIndiTrnovrPitGet(ticker=z1,secID="",endDate=tt,beginDate=t0,
                                           beginYear=u"",endYear=u"",reportType=u"",publishDateEnd=u"",
                                           publishDateBegin=u"",field=u"",pandas="1")
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn2_d1))
        #S19扩展
        """
        #7 单季度财务指标 DataAPI.FdmtIndiQGet(ticker=u"688002",secID="",endDate="",beginDate="",beginYear=u"",endYear=u"",reportType=u"",field=u"",pandas="1") #数据库yq_FdmtIndiQGet
        def get_FdmtIndiQGet():
            fn_d1= 'FdmtIndiQGet_%s' % tt
            fn1_d1 = '%s.csv' % fn_d1
            fn2_d1 = '%s.zip' % fn_d1
            info = '财务指标-单季度财务指标'
            if fn2_d1 not in list_files():
                x=DataAPI.FdmtIndiQGet(ticker=z1,secID="",endDate="",beginDate=t0,
                                       beginYear=u"",endYear=u"",reportType=u"",field=u"",pandas="1")
                save_data_adair(fn_d1,x)
                print('%s已经更新到%s' % (info,fn2_d1))  
            else:
                print('%s已经存在，%s数据已经更新，未执行' % (fn2_d1,info)) 
                
        """
        #7 财务指标—每股 (Point in time)
        #DataAPI.FdmtIndiPSPitGet(ticker=u"688002",secID="",endDate="",beginDate="",beginYear=u"",endYear=u"",reportType=u"",publishDateEnd=u"",publishDateBegin=u"",field=u"",pandas="1")
        def get_FdmtIndiPSPitGet():
            fn_d1= 'FdmtIndiPSPitGet_%s' % tt
            fn2_d2 = '%s.csv' % fn_d1            
            if fn2_d2 in list_files():
                return
            fn2_d1 = fn_d1
            info = '财务指标-每股'
            x=DataAPI.FdmtIndiPSPitGet(ticker=z1,secID="",endDate=tt,beginDate=t0,
                                           beginYear=u"",endYear=u"",reportType=u"",publishDateEnd=u"",
                                           publishDateBegin=u"",field=u"",pandas="1")
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn2_d1))  
            
        #8 合并现金流量表 DataAPI.FdmtCFGet(ticker=u"688001",secID=u"",reportType=u"",endDate=u"",beginDate=u"",
        # publishDateEnd=u"",publishDateBegin=u"",endDateRep="",beginDateRep="",beginYear="",endYear="",fiscalPeriod="",field=u"",pandas="1")
        def get_FdmtCFGet():
            fn_d1= 'FdmtCFGet_%s' % tt
            fn2_d2 = '%s.csv' % fn_d1            
            if fn2_d2 in list_files():
                return
            fn2_d1 = fn_d1
            info = '合并现金流量表'
            x=DataAPI.FdmtCFGet(ticker=z1,secID="",endDate=tt,beginDate=t0,
                                           beginYear=u"",endYear=u"",reportType=u"",publishDateEnd=u"",
                                           publishDateBegin=u"",field=u"",pandas="1")
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn2_d1))  
        #9 财务指标—盈利能力 (Point in time)
        # DataAPI.FdmtIndiRtnPitGet(ticker=u"688002",secID="",endDate="",beginDate="",
        #  beginYear=u"",endYear=u"",reportType=u"",publishDateEnd=u"",publishDateBegin=u"",field=u"",pandas="1")
        def get_FdmtIndiRtnPitGet():
            fn_d1= 'FdmtIndiRtnPitGet_%s' % tt
            fn2_d2 = '%s.csv' % fn_d1            
            if fn2_d2 in list_files():
                return
            fn2_d1 = fn_d1
            info = '财务指标—盈利能力'
            x=DataAPI.FdmtIndiRtnPitGet(ticker=z1,secID="",endDate=tt,beginDate=t0,
                                           beginYear=u"",endYear=u"",reportType=u"",publishDateEnd=u"",
                                           publishDateBegin=u"",field=u"",pandas="1")
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn2_d1))  
            
        #10 合并利润表TTM
        def FdmtISTTMPITGet_adair():
            fn_d1= 'FdmtISTTMPITGet_%s' % tt
            fn2_d2 = '%s.csv' % fn_d1            
            if fn2_d2 in list_files():
                return
            fn2_d1 = fn_d1
            info = '合并利润表TTM'
            x=DataAPI.FdmtISTTMPITGet(ticker=z1,secID=u"",endDate=tt,beginDate=t0,publishDateEnd=u"",publishDateBegin=u"",field=u"",pandas="1")        
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn2_d1))      
        
        #11 合并现金流量表（TTM Point in time）
        def FdmtCFTTMPITGet_adair():
            fn_d1= 'FdmtCFTTMPITGet_%s' % tt
            fn2_d2 = '%s.csv' % fn_d1            
            if fn2_d2 in list_files():
                return
            fn2_d1 = fn_d1
            info = '合并现金流量表（TTM Point in time）'
            x=DataAPI.FdmtCFTTMPITGet(ticker=z1,secID=u"",endDate=tt,beginDate=t0,publishDateEnd=u"",publishDateBegin=u"",field=u"",pandas="1")        
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn2_d1))  
        
        #12 合并利润表（单季度 Point in time）    
        #DataAPI.FdmtISQPITGet(ticker=u"688002",secID="",endDate="",beginDate="",beginYear=u"",endYear=u"",reportType=u"",publishDateEnd="",publishDateBegin="",isNew="",isCalc="",field=u"",pandas="1")    
        def FdmtISQPITGet_adair():
            fn_d1= 'FdmtISQPITGet_%s' % tt
            fn2_d2 = '%s.csv' % fn_d1            
            if fn2_d2 in list_files():
                return
            fn2_d1 = fn_d1
            info = '合并利润表TTM'
            x=DataAPI.FdmtISQPITGet(ticker=z1,secID="",endDate=tt,beginDate=t0,beginYear=u"",endYear=u"",reportType=u"",
                                    publishDateEnd=u"",publishDateBegin=u"",isNew="",isCalc="",field=u"",pandas="1")       
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn2_d1))   
        
        #14 业绩预告 S49-p1
        def get_FdmtEfGet_S49():
            fn_d1= 'FdmtEfGet_S49_%s' % tt
            fn2_d2 = '%s.csv' % fn_d1            
            if fn2_d2 in list_files():
                return
            fn2_d1 = fn_d1
            info = '业绩快报'
            x=DataAPI.FdmtEfGet(secID=u"",ticker=z1,publishDateBegin=u"",endDate=tt,beginDate=t0,field=u"",pandas="1")
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn2_d1))  
            
        #16 S53更新
        def get_FdmtISQGet():
            fn_d1= 'FdmtISQGetS53_%s' % tt
            fn2_d2 = '%s.csv' % fn_d1            
            if fn2_d2 in list_files():
                return
            fn2_d1 =  fn_d1
            info = '合并利润表（单季度，根据所有会计期末最新披露数据计算）'
            x=DataAPI.FdmtISQGet(ticker=u"",secID=u"",endDate=tt,beginDate=t0,beginYear=u"",endYear=u"",reportType=u"",field=u"",pandas="1")       
            save_data_adair(fn_d1,x)
            print('%s已经更新到%s' % (info,fn2_d1))
        
        def do_update(ind):
            if ind==0:
                #1合并利润表
                get_FdmtISGet()
            elif ind==1:
                #2资产重组表
                get_EquRestructuringGet()
            elif ind==2:
                #3合并资产负债表
                get_FdmtBSGet()
            elif ind==3:
                #4主营业务构成
                get_FdmtMainOperNGet()
            elif ind==4:
                #5财务衍生数据
                get_FdmtDerPitGet()
            elif ind==5:
                #6财务指标 营运能力
                get_FdmtIndiTrnovrPitGet()
            elif ind==6:
                #7 财务衍生数据
                get_FdmtIndiPSPitGet()
            elif ind==7:
                #8 合并现金流量表
                get_FdmtCFGet()
            elif ind==8:
                #9 财务指标—盈利能力 (Point in time)
                get_FdmtIndiRtnPitGet()
            elif ind==9:
                #10 合并利润表TTM
                FdmtISTTMPITGet_adair()
            elif ind==10:
                #11 合并现金流量表（TTM Point in time）
                FdmtCFTTMPITGet_adair()
            elif ind==11:
                #12 合并利润表（单季度 Point in time）    
                FdmtISQPITGet_adair()
            elif ind==12:
                #13 业绩快报
                get_FdmtEeGet_S26()
            elif ind==13:
                #14 业绩预告
                get_FdmtEfGet_S49()        
            elif ind==14:
                #16 合并利润表（单季度，根据所有会计期末最新披露数据计算）
                get_FdmtISQGet()
        
        w_n = 14        
        ind_p=range(w_n+1)   
        if para_sel:
            #多进程
            pool = ThreadPool(processes=w_n+1)
            pool.map(do_update, ind_p)
            pool.close()
            pool.join()
        else:
            for i in ind_p:
                do_update(i)
        
        
def get_all_date(tt):
    return pd.date_range('2000-01-01',tt).astype(str).tolist()

def update_data_temp(imputdata):
    t0,tt=imputdata
    uq_m.get_rontie_20d(t0,tt,False)
#执行依次，至少使用api60次
if __name__ == '__main__':
    #获取初始时间处理   
    uq_m = uq_methods()    
    tt = datetime.datetime.strftime( datetime.datetime.today(),'%Y%m%d')
    t_1y = '%d%s' % (int(tt[:4])-1,tt[4:])
    t_hy = datetime.datetime.strftime( datetime.datetime.today()-datetime.timedelta(60),'%Y%m%d')
    t_20d=datetime.datetime.strftime( datetime.datetime.today()-datetime.timedelta(20),'%Y%m%d')
    t_10d=datetime.datetime.strftime( datetime.datetime.today()-datetime.timedelta(10),'%Y%m%d')
    #table_name = 'yq_MktStockFactorsOneDayProGet'.lower()
    #t0 = get_ini_data('ResConSecIncomeGetS18'.lower(),'repForeTime')
    t0='2001-01-01'
    #t01 = get_ini_data('yq_MktStockFactorsOneDayProGet'.lower(),'tradeDate')
    #t02= get_ini_data('MktCmeFutdGet_S50'.lower(),'tradeDate')
    #t0=min(t01,t02)
    #tref0=get_tradingdate_adair(tt)
    tref0 = get_all_date(tt)
    tref0=[i.replace('-','') for i in tref0 if i>=t0]
    t0_3d,tt_3d=tref_split(tref0,3)
    t0_10d,tt_10d=tref_split(tref0,10)
    t0_20d,tt_20d=tref_split(tref0,20)
    
    _,fn_exist=get_file_name(datadir,'.csv')
    pool = ThreadPool(processes=num_core)
    pool.map(update_data_temp, zip(t0_20d,tt_20d))
    pool.close()
    pool.join()
    
    obj_t.use()