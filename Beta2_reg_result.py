import pandas as pd
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import json

# %%
with open('./component/500.json', 'r') as f:
    com_500 = json.load(f)
stk_list = glob.glob('./ind500_month/*.csv')

# %% 获取数据
year_data = {}
for year in com_500.keys():
    if len(com_500[year]) != 0:
        year_dfs = []
        for stk in com_500[year]:
            path = './ind500_month\\%s.csv' % stk.split('.')[0]
            if path in stk_list:
                df = pd.read_csv('./ind500_month\\%s.csv' % stk.split('.')[0], index_col=0)
                df = df.loc[(df.index > int(year) * 100) & (df.index < (int(year) + 1) * 100)]
                df['month'] = df.index.map(str)
                # df.set_index('ts_code', append=True, inplace=True)

                year_dfs.append(df)
        year_df = pd.concat(year_dfs)
        year_df.index = range(year_df.shape[0])

        year_data[year] = year_df


# %%

def read_index():
    df = pd.read_csv('./index_data/hs_500.csv', index_col=0)

    def split_date(x):
        y, m, d = x.split('-')
        return y + m + d

    df['trade_date'] = df['tradeDate'].apply(split_date)
    df['month'] = df['trade_date'].apply(lambda x: x[:6])
    return df.groupby('month')['CHGPct'].apply(lambda x: (x + 1).prod() - 1)


# %%
index_ser = read_index()
dummy_ser = index_ser.shift(1).apply(lambda x: 1 if x < 0 else 0)
dummy_ser.name = 'dummy'


# %%
def quantile_beta_diff(x, quantile=10):
    quan10 = np.quantile(x, quantile / 100)
    quan90 = np.quantile(x, 1 - quantile / 100)
    high, low = x[x > quan90].mean(), x[x < quan10].mean()
    return high - low


# %%
beta_ser = pd.concat([year_data[year].groupby('month')['beta'].apply(lambda x: quantile_beta_diff(x, quantile=10))
                      for year in year_data.keys()])
plt.figure(figsize=(16, 8))
plt.subplot(121)
plt.title('beta dispersion')
beta_ser.plot()
plt.subplot(122)
plt.title('index')
(index_ser.loc[beta_ser.index] + 1).cumprod().plot()
plt.show()

# %%
import statsmodels.api as sm

result = sm.OLS(index_ser.loc[beta_ser.index], sm.add_constant(beta_ser))
result = result.fit()

print(result.summary())

# %%
result = sm.OLS(index_ser.loc[beta_ser.index],
                pd.concat([beta_ser, dummy_ser.loc[beta_ser.index] * beta_ser], axis=1))
result = result.fit()

print(result.summary())

# %% strategy- beta > quantile and dummy = 0
signal = pd.concat([beta_ser, dummy_ser.loc[beta_ser.index]], axis=1)
quan_signal = pd.Series([beta_ser.iloc[:i].quantile() for i in range(beta_ser.shape[0] - 1)], index=signal.index[1:])
quan_signal.loc[signal.index[0]] = np.nan
quan_signal.sort_index(inplace= True)
# 当前时间段的beta阈值是从开始时间到上一个时间点时间段的beta序列的中位数。

signal = (signal['beta'] > quan_signal) & (signal['dummy'] == 0)
(index_ser.loc[beta_ser.index] + 1).cumprod().plot()
(index_ser.loc[beta_ser.index] * signal + 1).cumprod().plot()
plt.show()

# %% rolling predict
data_x = sm.add_constant(pd.concat([beta_ser, dummy_ser.loc[beta_ser.index]], axis=1))
data_y = index_ser.loc[beta_ser.index]

# %%
dfs = []
window = 35
for i in range(data_y.shape[0] - window):
    x = data_x.iloc[:i + window]
    y = data_y.iloc[:i + window]
    x_pre = data_x.iloc[[i + window]]
    model = sm.OLS(y, x).fit()
    y_pre = model.predict(x_pre)
    dfs.append(y_pre)
signal = pd.concat(dfs)
(index_ser.loc[signal.index] + 1).cumprod().plot()
((signal > 0) * index_ser.loc[signal.index] + 1).cumprod().plot()
plt.show()

# %%
