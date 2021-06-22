import pandas as pd
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from glob import glob

# %%
month_list = glob('./month_data/*.csv')
month_list.sort()

deep_list = glob('./deep_result/*.csv')
deep_list.sort()

dfs = {}
for path_ic in tqdm(month_list):
    date = path_ic[-11:-4]
    df = pd.read_csv(path_ic, index_col=0)
    df_ = df.iloc[:, 3:-3]
    ic_df = df_.corrwith(df_['return']).iloc[:-1]
    dfs[date] = (df, ic_df)

# %%
rets = {}
month_step = 6
for i in range(len(month_list) - month_step):
    train_dates = month_list[i: i + month_step]
    ics = []
    for path in train_dates:
        date = path[-11:-4]
        df, ic = dfs[date]
        ics.append(ic)
    ics = pd.concat(ics, axis=1).mean(axis=1)
    max_id = ics.abs().argmax()
    fac_id, sig = ics.index[max_id], np.sign(ics.iloc[max_id])
    test_date = month_list[i + month_step][-11:-4]
    test_df = dfs[test_date][0]
    deep_df = pd.read_csv(deep_list[i], index_col=0)
    deep_df['market_value'] = test_df['market_value']
    assert deep_list[i][-11:-4] == test_date

    if sig > 0:
        deep_df['max_ic'] = pd.qcut(test_df[fac_id], 3, [1, 2, 3])
    else:
        deep_df['max_ic'] = pd.qcut(test_df[fac_id], 3, [3, 2, 1])

    long_df = deep_df[deep_df['max_ic'] == 3].sort_values('pred', ascending=False).iloc[:35]
    short_df = deep_df[deep_df['max_ic'] == 1].sort_values('pred', ascending=True).iloc[:35]
    long_deep_df = deep_df.sort_values('pred', ascending=False).iloc[:35]
    short_deep_df = deep_df.sort_values('pred', ascending=True).iloc[:35]
    rets[test_date] = (long_df['real_return'].mean(),
                       short_df['real_return'].mean(),
                       long_deep_df['real_return'].mean(),
                       short_deep_df['real_return'].mean(),
                       )
    # rets[test_date] = ((long_df['real_return'] * (long_df['market_value'] / long_df['market_value'].sum())).sum(),
    #                    (short_df['real_return'] * (short_df['market_value'] / short_df['market_value'].sum())).sum(),
    #                    (long_deep_df['real_return'] * (long_deep_df['market_value'] / long_deep_df['market_value'].sum())).sum(),
    #                    (short_deep_df['real_return'] * (short_deep_df['market_value'] / short_deep_df['market_value'].sum())).sum(),
    #                    )

# %%
ret_df = pd.DataFrame(rets, index=['long', 'short', 'long_deep', 'short_deep']).T
(ret_df + 1).cumprod().plot(figsize=(16, 8))
plt.show()

# %%
(ret_df['long'] - ret_df['short'] + 1 - 0.005).cumprod().plot(figsize=(16, 8), label='enhence')
(ret_df['long_deep'] - ret_df['short_deep'] + 1 - 0.005).cumprod().plot(figsize=(16, 8), label='deep')
plt.legend()
plt.show()


# %% charge index 800
hs_800 = pd.read_csv('hsm_800.csv', index_col = 0)
hs_800['month'] = hs_800['endDate'].apply(lambda x: x[:7])
hs_800.index = hs_800['month']
ret_df['index'] = hs_800['chgPct'].shift(1)

# %%
(ret_df + 1).cumprod().plot(figsize=(16, 8))
plt.show()

# %%
(ret_df['long'] - ret_df['index'] + 1).cumprod().plot(figsize=(16, 8), label='enhence')
(ret_df['long_deep'] - ret_df['index'] + 1).cumprod().plot(figsize=(16, 8), label='deep')
plt.legend()
plt.show()