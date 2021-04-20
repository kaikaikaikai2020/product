import pandas as pd
import matplotlib.pyplot as plt
# %%
fac_df = pd.read_csv('fac_df.csv')

# %%
fac_df['qcut'] = fac_df.groupby('code')['1d_20d'].apply(lambda x: pd.qcut(x, 10, range(1, 11)))
fac_df['qcut_5'] = fac_df.groupby('code')['5d_20d'].apply(lambda x: pd.qcut(x, 10, range(1, 11)))
fac_df['qcut_20'] = fac_df.groupby('code')['20d_20d'].apply(lambda x: pd.qcut(x, 10, range(1, 11)))

# %%
plt.figure(figsize=(16,8))
ret = fac_df.groupby(['date', 'qcut'])['chg_pct'].mean().unstack()
(ret[1] - ret[10] + 1).cumprod().plot(label = '1d_20d')

ret = fac_df.groupby(['date', 'qcut_5'])['chg_pct'].mean().unstack()
(ret[1] - ret[10] + 1).cumprod().plot(label = '5d_20d')

ret = fac_df.groupby(['date', 'qcut_20'])['chg_pct'].mean().unstack()
(ret[1] - ret[10] + 1).cumprod().plot(label = '20d_20d')
plt.legend()
plt.show()