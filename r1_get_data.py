import pandas as pd
import numpy as np
from glob import glob
import os
from tqdm import tqdm
import json

# %%
stk_list = glob('./new_stk_data/new_stk_data/*/*.*.CSV')
use_cols = ['代码', '简称', '日期', '开盘价(元)', '最高价(元)',
            '最低价(元)', '收盘价(元)', '成交量(股)',
            '换手率(%)', 'A股流通市值(元)', '总市值(元)', '市盈率']
cols_name = ['ts_code', 'short_name', 'trade_date', 'open', 'high',
             'low', 'close', 'vol',
             'turnover', 'flowvol', 'value', 'pb']


# %%
def read_file(path):
    df = pd.read_csv(path, encoding='gbk', usecols=use_cols, na_values='--')
    df.columns = cols_name
    df.dropna(inplace=True)
    flag = True
    if df.shape[0] < 200:
        flag = False

    def split_date(x):
        y, m, d = x.split('-')
        return y + m + d

    if flag:
        df['trade_date'] = df['trade_date'].apply(split_date)
        df['pre_close'] = df['close'].shift(1)
        df['chg_pct'] = df['close'] / df['pre_close'] - 1
        df.sort_values('trade_date', inplace=True)
        return df.iloc[60:], flag
    else:
        return None, flag


def read_processing_stk_data(path):
    df, flag = read_file(path)
    if flag:
        df['month'] = df['trade_date'].apply(lambda x: x[:6])
        df.index = df.pop('trade_date')
        df['ts_code'] = df['ts_code'].apply(lambda x: x.split('.')[0])
        df.sort_index(inplace=True)
        df_month = df.groupby('month').agg({'chg_pct': lambda x: (x + 1).prod() - 1,
                                            'close': lambda x: x.iloc[-1],
                                            'open': lambda x: x.iloc[0],
                                            'pb': lambda x: x.mean(),
                                            'value': lambda x: x.mean()})
        ts_code = df['ts_code'].unique()[0]
        return df_month, ts_code
    else:
        return None, flag


# %%
df, _ = read_processing_stk_data(stk_list[1])

# %%
import multiprocessing

if __name__ == '__main__':
    num_core = 8
    file_step = 16
    n = int(len(stk_list) / 16 + 1)

    for i in tqdm(range(n)):
        with multiprocessing.Pool(num_core) as pool:
            s = pool.map(read_processing_stk_data, stk_list[i * file_step: (i + 1) * file_step])
            for (ret_df, ts_code) in s:
                if ts_code:
                    if not os.path.exists('./month_data/'): os.mkdir('./month_data')
                    ret_df.to_csv('./month_data/%s.csv' % ts_code)
