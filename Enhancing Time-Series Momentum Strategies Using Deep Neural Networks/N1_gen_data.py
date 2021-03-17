import os
import pandas as pd
import numpy as np

# %% Features
def rocp_func(c):
    c_pre = c.shift(1)
    return (c - c_pre) / c_pre


def orocp_func(o):
    o_pre = o.shift(1)
    return (o - o_pre) / o_pre


def hrocp_func(h):
    h_pre = h.shift(1)
    return (h - h_pre) / h_pre


def lrocp_func(l):
    l_pre = l.shift(1)
    return (l - l_pre) / l_pre


def ma_func(c, n=5):
    return c.rolling(n).mean()


def macp_func(ma, c):
    return (ma - c) / c


def ewma_func(c, n=5):  # 传统ewma计算方法会导致时间长度不一致产生的数值不同问题，按照加权逻辑自定义一种计算过程
    base = (n - 1) / (n + 1)
    w_ = np.array([base ** i for i in range(1, n + 1)])
    w_ = w_ / w_.sum()
    return c.rolling(n).apply(lambda x: (x * w_).sum(), raw=True)


def macd_func(c, n_fast=5, n_slow=12, n_signal=6):
    ewmafast = ewma_func(c, n=n_fast)
    ewmaslow = ewma_func(c, n=n_slow)
    dif = ewmafast - ewmaslow
    dea = ewma_func(dif, n=n_signal)
    macd_ = (dif - dea)
    return pd.DataFrame([dif, dea, macd_], index=['dif', 'dea', 'macd'], columns=c.index).T


def difrocp_func(dif):
    dif_pre = dif.shift(1)
    return (dif - dif_pre) / dif_pre


def dearocp_func(dea):
    dea_pre = dea.shift(1)
    return (dea - dea_pre) / dea_pre


def macdrocp_func(macd):
    macd_pre = macd.shift(1)
    return (macd - macd_pre) / macd_pre


# def rsi_func(c, n=12): # 出现前后不一致
#     c_v = c.values
#     rsi = talib.RSI(c_v, timeperiod=n)
#     return pd.Series(rsi, index=c.index)

def rsi_func(c, n=12):
    r = c / c.shift(1) - 1
    up_rate = r.rolling(n).apply(lambda x: x[x > 0].sum() / x.shape[0], raw=True)
    down_rate = r.rolling(n).apply(lambda x: -x[x < 0].sum() / x.shape[0], raw=True)
    rs = up_rate / down_rate
    # 异常值处理
    rsi = 100 - 100 / (1 + rs)
    rsi.replace(np.inf, 1, inplace=True)
    rsi.replace(0, 0.1, inplace=True)
    return rsi


def rsirocp_func(rsi):
    rsi_pre = rsi.shift(1)
    return (rsi - rsi_pre) / rsi_pre


def vrocp_func(v):
    v_pre = v.shift(1)
    return np.arctan((v - v_pre) / v_pre)

def boll_fn(cn, timeperiod = 20):
    ma = ma_func(cn, n = timeperiod)
    std_ = cn.rolling(timeperiod).std()
    return ma + 2 * std_, ma, ma - 2 * std_
    

def boll_func(c, n=6):
    # c_v = c.copy()
    upper, middle, lower = boll_fn(c, n)
    return pd.DataFrame([upper, middle, lower], index=['up', 'mid', 'low'], columns=c.index).T


def vmarocp_func(v, n=5):
    v_ma = v.rolling(n).mean()
    v_pre = v_ma.shift(1)
    return (v_ma - v_pre) / v_pre


def vmacp_func(v, n=5):
    v_ma = v.rolling(n).mean()
    return (v_ma - v) / v


def rocpvrocp_func(rocp, vrocp):
    return rocp * vrocp


def main_features(df_ori):
    c, l, o, h, v = [df_ori[col] for col in ['closePrice', 'lowestPrice', 'openPrice', 'highestPrice', 'turnoverVol']]
    rocp = rocp_func(c)
    orocp = orocp_func(o)
    lrocp = lrocp_func(l)
    hrocp = hrocp_func(h)
    ma = ma_func(c)
    macp = macp_func(ma, c)
    macd_df = macd_func(c)
    dif, dea, macd = [macd_df[col] for col in macd_df.columns]
    difrocp = difrocp_func(dif)
    dearocp = dearocp_func(dea)
    macdrocp = macdrocp_func(macd)
    rsi = rsi_func(c)
    rsirocp = rsirocp_func(rsi)
    vrocp = vrocp_func(v)
    boll_df = boll_func(c)
    up, mid, low = [boll_df[col] for col in boll_df.columns]
    vmarocp = vmarocp_func(v)
    vmacp = vmacp_func(v)
    rocpvrocp = rocpvrocp_func(rocp, vrocp)

    df_feature = pd.concat([rocp, orocp, lrocp, hrocp, ma, macp, dif, dea, macd, difrocp, dearocp,
                            macdrocp, rsi, rsirocp, vrocp, up, mid, low, vmarocp, vmacp, rocpvrocp], axis=1)
    cols_name = ['rocp', 'orocp', 'lrocp', 'hrocp', 'ma', 'macp', 'dif', 'dea', 'macd', 'difrocp', 'dearocp',
                 'macdrocp', 'rsi', 'rsirocp', 'vrocp', 'up', 'mid', 'low', 'vmarocp', 'vmacp', 'rocpvrocp']
    df_feature.columns = cols_name
    # 添加y的信息
    df_feature.index = df_ori.index
    df_feature['sigma'] = df_ori['chgPct'].rolling(252).std()
    df_feature['code'] = df_ori['contractObject']
    df_feature['chgPct'] = df_ori['chgPct']
    df_feature['nextPct'] = df_ori['chgPct'].shift(-1)
    df_feature['preClosePrice'] = df_ori['preClosePrice']
    df_feature['closePrice'] = df_ori['closePrice']
    df_feature['nextPrice'] = df_ori['nextPrice']
    df_feature.index = df_ori['tradeDate']
    return df_feature


def read_ticker_name(name='SH600000'):
    df_ori = pd.read_csv('./stk_data/%s.csv' % name)
    df_ori['nextPrice'] = df_ori['closePrice'].shift(-1)
    df_ori.index = np.arange(df_ori.shape[0])
    return df_ori


def get_data_by_date(date, data_list):
    def get_data(data, date):
        if date in data.index:
            return data.loc[date]

    return [get_data(data, date) for data in data_list]


# %% 读取数据
df_ori = pd.read_csv('./data/yq_MktMFutdGet.csv')
tickers = set(df_ori['contractObject'].tolist())
dates = set(df_ori['tradeDate'].tolist())

# %% 分割数据
for sym in tickers:
    sub_df = df_ori[df_ori['contractObject'] == sym]
    sub_df.sort_values('tradeDate', inplace=True)
    if not os.path.exists('./stk_data'): os.mkdir('./stk_data')
    sub_df.to_csv('./stk_data/%s.csv'%sym, index=False)

# %% 计算特征
feat_list = []
for ind in tickers:
    df_ticker = pd.read_csv('./stk_data/%s.csv'%ind)
    df_ticker.index = np.arange(df_ticker.shape[0])
    df_ticker['nextPrice'] = df_ticker['closePrice'].shift(-1)
    df_feature = main_features(df_ticker)
    feat_list.append(df_feature.dropna())
    print(ind)
    if not os.path.exists('./feat_data'): os.mkdir('./feat_data')
    df_feature.to_csv('./feat_data/%s.csv' % ind)

# %% 整合成时间文件
df_time = pd.concat(feat_list)
for ind in set(df_time.index.tolist()):
    if not os.path.exists('./date_data'): os.mkdir('./date_data')
    df_time.loc[ind].to_csv('./date_data/%s.csv' % ind)
