import pandas as pd
import numpy as np
import os
from glob import glob
from tqdm import tqdm


# %%
def read_file(path):
    df = pd.read_csv(path, dtype={'month': str})
    return df[(df['month'] > '2000') & (df['month'] < '2021')]


def read_index():
    df = pd.read_csv('./index_data/sh_index.csv', index_col=0)

    def split_date(x):
        y, m, d = x.split('-')
        return y + m + d

    df['trade_date'] = df['tradeDate'].apply(split_date)
    df['month'] = df['trade_date'].apply(lambda x: x[:6])
    return df.groupby('month')[['CHGPct']].apply(lambda x: (x + 1).prod() - 1)


# %%
stk_list = glob('./month_data/*.csv')
data = pd.concat([read_file(path) for path in stk_list])
df_index = read_index()

# %%
month_list = data['month'].unique().tolist()
month_list.sort()

# %%
fac_dict = {}
for month in tqdm(month_list):
    month_data = data[data['month'] == month].copy()
    sb = pd.qcut(month_data['value'], 2, ['s', 'b'])
    hml = pd.qcut(month_data['pb'], 3, ['l', 'm', 'h'])

    fac_df = month_data.groupby([sb, hml]).apply(
        lambda x: (x['value'] / x['value'].sum() * x['chg_pct']).sum()).unstack()

    smb_ = (fac_df.loc['s'] - fac_df.loc['b']).mean()
    hml_ = (fac_df.loc[:, 'l'] - fac_df.loc[:, 'h']).mean()
    fac_dict[month] = {'smb': smb_, 'hml': hml_}

# %%
fac_df = pd.concat([df_index, pd.DataFrame(fac_dict).T], axis=1).dropna()

# %%
for path in tqdm(stk_list):
    df = read_file(path)
    df.set_index('month', inplace=True)
    df = pd.concat([df, fac_df], axis=1).dropna()
    if df.shape[0] > 36:
        if not os.path.exists('./fac_data/'): os.mkdir('./fac_data/')
        df.to_csv('./fac_data/%s' % path.split('\\')[-1])
