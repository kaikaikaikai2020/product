import pandas as pd
import numpy as np
from glob import glob
import os
from tqdm import tqdm

# %%
base_dir = '/home/dinker/桌面/codes/DayNightChina/new_stk_data'
stk_list = glob(base_dir + '/*/*.CSV')

# %%
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
        df['month'] = df['trade_date'].apply(lambda x: x[:6])
        df.sort_values('trade_date', inplace=True)
        return df.iloc[60:], flag
    else:
        return None, flag


# %%
# df = read_file(stk_list[1])[0]


# %%
# ind = df['month'].tolist()
#
# # %%
# a = df['chg_pct']
#
# a1, a2, a3 = a[:40], a[40:80], a[:80]
# print(a1.std(), a2.std(), a3.std(),
#       np.sqrt((a1.var() + a2.var()) / 2))

# %%
# df_month_ret = df.groupby('month')['chg_pct'].apply(lambda x: (1 + x).prod() - 1)
# df_month_ret.rolling(3).apply(lambda x: (1 + x).prod() - 1)


# %% 计算macd
def macd_func(df, s, l, p_window=63, r_window=252):
    df = df.copy()
    y_ks = []
    w = 1 / len(s)
    y_s = pd.Series(0, index=df.index)
    for s_k, l_k in zip(s, l):
        short = df['close'].ewm(adjust=False, halflife=np.log(0.5) / np.log(1 - 1 / s_k), ignore_na=True).mean()
        long = df['close'].ewm(adjust=False, halflife=np.log(0.5) / np.log(1 - 1 / l_k), ignore_na=True).mean()
        macd = short - long
        eta = macd / df['close'].rolling(p_window).std()  # y_k
        y_k = eta / eta.rolling(r_window).std()  # standar y_k
        y_ks.append(y_k)
        y_s += w * y_k * y_k.apply(lambda x: -x ** 2 / 4) / 0.89
        df['macd_%s_%s' % (s_k, l_k)] = y_k
    df['macd_mom'] = y_s
    return df


def get_16_pred_sig(macd_df, windows=[2, 6, 12]):
    macd_cols = [i for i in macd_df.columns if 'macd' in i]
    base_macd = macd_df.groupby('month')[macd_cols].mean()
    base_macd.columns = [col + '_1' for col in base_macd.columns]
    group_datas = [base_macd]
    for i in windows:
        g_df = base_macd.rolling(i).mean()
        g_df.columns = [col + '_%s' % i for col in g_df.columns]
        group_datas.append(g_df)
    return pd.concat(group_datas, axis=1)


def get_y(df):
    df_month_ret = df.groupby('month')['chg_pct'].apply(lambda x: (1 + x).prod() - 1)
    df_month_ret.name = 'next_month_ret'
    return df_month_ret.shift(-1)
    # df_month_short = df.groupby('month')['chg_pct'].apply(lambda x: (1 - x).prod() - 1)
    # df_month_short.name = 'next_month_short'
    # return pd.concat([df_month_ret.shift(-1), df_month_short.shift(-1)], axis=1)


def main_func(path):
    df, flag = read_file(path)
    # print(flag)
    if flag:
        # print(flag)
        ts_code = df['ts_code'].unique()[0].split('.')[0]
        macd_df = macd_func(df, [8, 16, 32], [24, 48, 96])
        if macd_df.shape[0] < 12:  # make sure has data
            ts_code = False
        x_df = get_16_pred_sig(macd_df, windows=[2, 6, 12])
        y_df = get_y(df)
        return pd.concat([x_df, y_df], axis=1).dropna(), ts_code
    else:
        return None, flag


def read_stk(path):
    ts_code = path[-10:-4]
    df = pd.read_csv(path, dtype={'month': str})
    df['ts_code'] = ts_code
    return df


# %%
import multiprocessing

# for stk_path in stk_list:
#     df, ts_code = main_func(stk_path)
#     # print(df)
#     if ts_code:
#         if not os.path.exists('./month_data/'): os.mkdir('./month_data')
#         df.to_csv('./month_data/%s.csv' % ts_code)

if __name__ == '__main__':
    num_core = 8
    file_step = 16
    n = int(len(stk_list) / 16 + 1)

    # saving month factor data
    for i in tqdm(range(n)):
        with multiprocessing.Pool(num_core) as pool:
            s = pool.map(main_func, stk_list[i * file_step: (i + 1) * file_step])
            for (ret_df, ts_code) in s:
                if ts_code:
                    if not os.path.exists('./month_data/'): os.mkdir('./month_data')
                    ret_df.to_csv('./month_data/%s.csv' % ts_code)

    # saving month date data
    stk_list = glob('./month_data/*.csv')

    df = pd.concat([read_stk(path) for path in stk_list])
    df.set_index(['month', 'ts_code'], inplace=True)
    month_list = df.index.get_level_values('month').unique().sort_values().to_numpy()
    month_list = month_list[month_list > '2005']
    if not os.path.exists('./month'): os.mkdir('./month')

    for month in tqdm(month_list):
        sub_df = df.loc[month]
        sub_df = sub_df.sort_values('next_month_ret', ascending=False)
        sub_df.to_csv('./month/%s.csv' % month)

# %%
# df = read_file(stk_list[3])
# macd_df = macd_func(df, [8, 16, 32], [24, 48, 96])
# x_df = get_16_pred_sig(macd_df, windows=[2, 6, 12])
# y_df = get_y(df)

# %%
