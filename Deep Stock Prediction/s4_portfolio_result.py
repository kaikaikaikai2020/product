import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

# %%
data_list = os.listdir('./pred_data')
stock_list = []
# %%
for a_type in ['lstm', 'lstmcm', 'arima']:
    ret_dfs = []
    sig_dfs = []
    for stk in ['GM', 'F', 'TM', 'TSLA']:
        df = pd.read_csv('./pred_data/%s_%s.csv' % (a_type ,stk), index_col=0)
        df['signal_%s' % stk] = df['pred'] / df['PRE_CLOSE'] - 1
        df['ret_%s' % stk] = df['CLOSE'] / df['PRE_CLOSE'] - 1
        ret_dfs.append(df['ret_%s' % stk])
        sig_dfs.append(df['signal_%s' % stk])
    
    # portfolio plot
    sig_df = pd.concat(sig_dfs, axis=1)
    ret_df = pd.concat(ret_dfs, axis=1)
    bool1 = (sig_df.T >= sig_df.max(axis=1)).T.values
    bool2 = (sig_df >= 0).values
    ret_arr = (ret_df.values * bool1 * bool2).sum(axis=1)
    cum_arr = np.cumprod(ret_arr + 1)
    plt.plot(cum_arr, label=a_type)
plt.title('with TSLA')
plt.legend()
plt.show()

# %%
for a_type in ['lstm', 'lstmcm', 'arima']:
    ret_dfs = []
    sig_dfs = []
    for stk in ['GM', 'F', 'TM']:
        df = pd.read_csv('./pred_data/%s_%s.csv' % (a_type, stk), index_col=0)
        df['signal_%s' % stk] = df['pred'] / df['PRE_CLOSE'] - 1
        df['ret_%s' % stk] = df['CLOSE'] / df['PRE_CLOSE'] - 1
        ret_dfs.append(df['ret_%s' % stk])
        sig_dfs.append(df['signal_%s' % stk])

    # portfolio plot
    sig_df = pd.concat(sig_dfs, axis=1)
    ret_df = pd.concat(ret_dfs, axis=1)
    bool1 = (sig_df.T >= sig_df.max(axis=1)).T.values
    bool2 = (sig_df >= 0).values
    ret_arr = (ret_df.values * bool1 * bool2).sum(axis=1)
    cum_arr = np.cumprod(ret_arr + 1)
    plt.plot(cum_arr, label=a_type)
plt.title('without TSLA')
plt.legend()
plt.show()