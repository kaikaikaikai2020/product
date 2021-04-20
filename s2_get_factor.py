import pandas as pd
import numpy as np
from glob import glob
from tqdm import tqdm
import matplotlib.pyplot as plt

# %%
ms = glob('./date_data/*.csv')


def read_file(path):
    df = pd.read_csv(path)
    df['date'] = path.split('.')[1][-10:]
    return df


# %%
datas = [read_file(path) for path in ms]

# %% get property data
df = pd.concat(datas)
df.set_index(['code', 'date'], inplace=True)
df.sort_index(level='code', inplace=True)

# %% split data to stock way
codes = df.index.get_level_values('code').unique().tolist()


# %%
def get_factor_data(df):
    sub_df = df.copy()
    sub_up = sub_df['high'] / sub_df['close'].shift(1) - 1
    sub_down = sub_df['low'] / sub_df['close'].shift(1) - 1
    # 将可能触及涨跌停的点去掉
    sub_df = sub_df[(sub_up < 0.1) & (sub_down > -0.1)]
    # 计算因子
    sub_df['chg_pct'] = sub_df['close'] / sub_df['close'].shift(1) - 1
    sub_df.dropna(inplace=True)
    sub_df['1d_20d'] = sub_df['arpp'].rolling(20).mean().shift(1)
    # 5d_20d
    high_5 = sub_df['high'].rolling(5).mean()
    low_5 = sub_df['low'].rolling(5).mean()
    appr_5d = sub_df['twap'].rolling(5).mean()
    sub_df['appr_5d'] = (appr_5d - low_5) / (high_5 - low_5)
    sub_df['5d_20d'] = sub_df['appr_5d'].rolling(20).mean().shift(1)
    # 20d_20d
    high_20 = sub_df['high'].rolling(20).mean()
    low_20 = sub_df['low'].rolling(20).mean()
    appr_20d = sub_df['twap'].rolling(20).mean()
    sub_df['appr_20d'] = (appr_20d - low_20) / (high_20 - low_20)
    sub_df['20d_20d'] = sub_df['appr_20d'].rolling(20).mean().shift(1)
    return sub_df.dropna()


# %% processing data
if __name__ == '__main__':
    factor_data = []
    for code in tqdm(codes):
        sub_df = df.loc[code].copy()
        sub_df['code'] = code
        sub_df.sort_index(inplace=True)
        if sub_df.shape[0] > 100:  # 保证数据的数量
            factor_data.append(sub_df)
    
    from multiprocessing import Pool
    
    fac_dfs = []
    n_step = 16
    n_round = int(len(ms) / n_step) + 1
    for i in tqdm(range(n_round)):
        with Pool(8) as pool:
            s = pool.map(get_factor_data, factor_data[i * n_step: (i + 1) * n_step])
            for df in s:
                fac_dfs.append(df)

    fac_df = pd.concat(fac_dfs)
    fac_df['date'] = fac_df.index
    fac_df.to_csv('fac_df.csv')


# %%
# fac_df = pd.concat(fac_dfs)
# fac_df.to_csv('fac_df.csv')
# fac_df['date'] = fac_df.index
# fac_df.index = range(fac_df.shape[0])
