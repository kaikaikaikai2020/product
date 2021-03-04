#日度框架测试
#HA配对交易
#tF '1.4.0'
#补充测试数据
import pandas as pd
import numpy as np
import MAIN.Basics as basics
import MAIN.Reinforcement as RL
import tensorflow as tf

#import seaborn as sns
#import matplotlib.pyplot as plt
from UTIL import FileIO
from STRATEGY.Cointegration import EGCointegration
from STRATEGY.Cointegration_pred import EGCointegration as EGCointegration_pred

from yq_toolsS45 import save_pickle,read_pickle
from yq_toolsS45 import list_intersec
from yq_toolsS45 import time_use_tool
import os
import multiprocessing
num_core = multiprocessing.cpu_count()
num_core = int(num_core/10*9)
pn0 = r'S60P3_para_CSI300'
if not os.path.exists(pn0):
    os.mkdir(pn0)
# Read config
config_path  = 'CONFIG\config_train.yml'
config_train = FileIO.read_yaml(config_path)


def updateS60data():
    obj_clock=time_use_tool()
    def trans_sub_x(sub_x):
        sub_x = sub_x.copy()
        sub_x.dropna(how='all',inplace=True)
        sub_x.dropna(how='all',axis='columns',inplace=True)
        sub_x=sub_x.iloc[2:]
        sub_x.columns=['date','close']
        sub_x['date'] = sub_x.date.apply(lambda x:x.date()).astype(str)
        return sub_x
    
    x = pd.read_excel('pair trading data CN HK KR TW.xlsx')
    x.dropna(how='all',inplace=True)
    x.dropna(how='all',axis='columns',inplace=True)
    x=x.iloc[2:]
    x.iloc[0].fillna(method='bfill',inplace=True)
    x.reset_index(drop =True,inplace=True)
    x.drop(1,inplace=True)
    x=x.T
     
    data = {}   
    while len(x)>0:
        sub_x = x.iloc[0:4].T
        x=x.iloc[4:]
        
        sub_x1=sub_x[sub_x.columns[0:2]]
        sub_x1 = trans_sub_x(sub_x1)
        sub_x2=sub_x[sub_x.columns[2:4]]
        sub_x2 = trans_sub_x(sub_x2)
        sub_ticker=sub_x.iloc[0].unique().tolist()
        sub_ticker = [i.split(' ')[0] for i in sub_ticker]
        data['-'.join(sub_ticker)] = [sub_x1,sub_x2]
    obj_clock.use('Time used')
    return data

def cal_HA_data(fn1,fn2,x=[],y=[],train_len=None):
    if train_len is None:
        com_fn1 =os.path.join(pn0, 'bac-%s-%s.csv' % (fn1,fn2))
        com_fn2 =os.path.join(pn0, 'par-%s-%s.pkl' % (fn1,fn2))
        com_fn3 =os.path.join(pn0, 'position-%s-%s.xlsx' % (fn1,fn2))
    else:
        com_fn1 =os.path.join(pn0, 'bac-%s-%s-%s.csv' % (fn1,fn2,train_len))
        com_fn2 =os.path.join(pn0, 'par-%s-%s-%s.pkl' % (fn1,fn2,train_len))
        com_fn3 =os.path.join(pn0, 'position-%s-%s-%s.xlsx' % (fn1,fn2,train_len))
    # Read prices
    if len(x)==0:
        x = pd.read_csv(r'STATICS\PRICE\%s.csv' % fn1)
        y = pd.read_csv(r'STATICS\PRICE\%s.csv' % fn2)
    tmp = list_intersec(x.date.tolist(),y.date.tolist())
    x=x[x.date.isin(tmp)]
    y=y[y.date.isin(tmp)]
    
    x, y = EGCointegration.clean_data(x, y, 'date', 'close')
    
    # Separate training and testing sets
    if train_len is None:
        train_pct = 0.5
        train_len = round(len(x) * train_pct)
    else:
        if not isinstance(train_len,int):
            train_len=x.date.tolist().index(train_len)
    idx_train = list(range(0, train_len))
    idx_test  = list(range(0, len(x)))
    EG_Train = EGCointegration(x.iloc[idx_train, :], y.iloc[idx_train, :], 'date', 'close')
    EG_Test  = EGCointegration_pred(x.iloc[idx_test,  :], y.iloc[idx_test,  :], 'date', 'close')
    
    # Create action space
    n_hist    = list(np.arange(20, 201, 20))
    n_forward = list(np.arange(20, 201, 20))
    trade_th  = list(np.arange(1,  5.1, 1))
    stop_loss = list(np.arange(1,  2.1, 0.5))
    cl        = list(np.arange(0.05,  0.11, 0.05))
    actions   = {'n_hist':    n_hist,
                 'n_forward': n_forward,
                 'trade_th':  trade_th,
                 'stop_loss': stop_loss,
                 'cl':        cl}
    n_action  = int(np.product([len(actions[key]) for key in actions.keys()]))
    
    # Create state space
    transaction_cost = [0.001]
    states  = {'transaction_cost': transaction_cost}
    n_state = len(states)
    
    # Assign state and action spaces to config
    config_train['StateSpaceState'] = states
    config_train['ActionSpaceAction'] = actions
    
    # Create and build network
    one_hot  = {'one_hot': {'func_name':  'one_hot',
                            'input_arg':  'indices',
                             'layer_para': {'indices': None,
                                            'depth': n_state}}}
    output_layer = {'final': {'func_name':  'fully_connected',
                              'input_arg':  'inputs',
                              'layer_para': {'inputs': None,
                                             'num_outputs': n_action,
                                             'biases_initializer': None,
                                             'activation_fn': tf.nn.relu,
                                             'weights_initializer': tf.ones_initializer()}}}
    
    state_in = tf.placeholder(shape=[1], dtype=tf.int32)
    
    N = basics.Network(state_in)
    N.build_layers(one_hot)
    N.add_layer_duplicates(output_layer, 1)
    
    # Create learning object and perform training
    RL_Train = RL.ContextualBandit(N, config_train, EG_Train)
    
    sess = tf.Session()
    if not os.path.exists(com_fn2):
        RL_Train.process(sess, save=False, restore=False)
        
        # Extract training results
        action = RL_Train.recorder.record['NETWORK_ACTION']
        reward = RL_Train.recorder.record['ENGINE_REWARD']
        print(np.mean(reward))
        
        df1 = pd.DataFrame()
        df1['action'] = action
        df1['reward'] = reward
        #mean_reward = df1.groupby('action').mean()
        #sns.distplot(mean_reward)
        
        # Test by trading continuously
        [opt_action] = sess.run([RL_Train.output], feed_dict=RL_Train.feed_dict)
        opt_action = np.argmax(opt_action)
        action_dict = RL_Train.action_space.convert(opt_action, 'index_to_dict')
        save_pickle(com_fn2,action_dict)
    else:
        action_dict=read_pickle(com_fn2)
    indices = range(action_dict['n_hist'], len(EG_Test.x)-1)
    
    pnl = pd.DataFrame()
    pnl['Time'] = EG_Test.timestamp
    pnl['Trade_Profit'] = 0
    pnl['Cost'] = 0
    pnl['N_Trade'] = 0
    
    import warnings
    warnings.filterwarnings('ignore')
    rec_mark = []
    for i in indices:
        if i % 100 == 0:
            print(i)
        EG_Test.process(index=i, transaction_cost=0.001, **action_dict)
        trade_record = EG_Test.record
        if (trade_record is not None) and (len(trade_record) > 0):
            print('value at {}'.format(i))
            trade_record = pd.DataFrame(trade_record)
            trade_cost   = trade_record.groupby('trade_time')['trade_cost'].sum()
            close_cost   = trade_record.groupby('close_time')['close_cost'].sum()
            profit       = trade_record.groupby('close_time')['profit'].sum()
            open_pos     = trade_record.groupby('trade_time')['long_short'].sum()
            close_pos    = trade_record.groupby('close_time')['long_short'].sum() * -1
    
            pnl['Cost'].loc[pnl['Time'].isin(trade_cost.index)] += trade_cost.values
            pnl['Cost'].loc[pnl['Time'].isin(close_cost.index)] += close_cost.values
            pnl['Trade_Profit'].loc[pnl['Time'].isin(close_cost.index)] += profit.values
            pnl['N_Trade'].loc[pnl['Time'].isin(trade_cost.index)] += open_pos.values
            pnl['N_Trade'].loc[pnl['Time'].isin(close_cost.index)] += close_pos.values
            trade_record['signal_time']=EG_Test.timestamp[i]
            rec_mark.append(trade_record)
    rec_mark = pd.concat(rec_mark)
    warnings.filterwarnings(action='once')   
    sess.close()
    #save result
    pnl.to_csv(com_fn1)
    #rec_mark.drop_duplicates(subset=['trade_time','close_time'],keep='last',inplace=True)
    #rec_mark.sort_values(by=['trade_time','close_time'],axis=0,ascending=[True,True],inplace=True)
    rec_mark.to_excel(com_fn3)
    
def do_cal(inputdata):
    fn1,fn2,x,y,train_len,check_sel = inputdata
    if train_len is not None:
        com_fn2 =os.path.join(pn0, 'par-%s-%s-%s.pkl' % (fn1,fn2,train_len))
        com_fn3 =os.path.join(pn0, 'Errorpar-%s-%s-%s.pkl' % (fn1,fn2,train_len))
    else:            
        com_fn2 =os.path.join(pn0, 'par-%s-%s.pkl' % (fn1,fn2))
        com_fn3 =os.path.join(pn0, 'Errorpar-%s-%s.pkl' % (fn1,fn2))
    if (os.path.exists(com_fn2) or os.path.exists(com_fn3)) and check_sel:
        return
    try:
        cal_HA_data(fn1,fn2,x,y,train_len)
    except:
        print('Error %s-%s' % (fn1,fn2))
        save_pickle(com_fn3,[])

"""
if __name__ == '__main__':
    f1 =['NTN_1M','NTN_1M','NTN_1M','NKY','TAMSCI','2330','HSI','XIN9I']
    f2 =['KWN_1M','NTN_6M','NTN_12M','TPX','TWSE','005930','HSCEI','HSCEI']
    f3 = [False]*len(f1)
    fn_a = zip(f1,f2,f3)
    pool = multiprocessing.Pool(num_core) 
    pool.map(do_cal, fn_a)
    pool.close()
    pool.join() 
"""

