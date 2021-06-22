from glob import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# %%
month_list = glob('./month_data/*.csv')
month_list.sort()

# %%
# path = month_list[0]
# df_ = df.iloc[:, 3:-3]

# %%
def get_condi_sub(df):
    condi_roe = df['ROE'] > np.quantile(df['ROE'], 0.75)
    condi_gpm = df['GrossIncomeRatio'] > np.quantile(df['GrossIncomeRatio'], 0.5)
    condi_yor = df['OperatingRevenueGrowRate'] > np.quantile(df['OperatingRevenueGrowRate'], 0.75)
    condi_yoa = df['TotalAssetGrowRate'] > np.quantile(df['TotalAssetGrowRate'], 0.5)
    condi_cfp = df['PCF'] < np.quantile(df['PCF'], 0.7)
    print(df[condi_roe & condi_gpm & condi_yor & condi_yoa & condi_cfp].shape)
    return df[condi_roe & condi_gpm & condi_yor & condi_yoa & condi_cfp]


# %%
rets = {}
for path in month_list:
    date = path[-11:-4]
    df = pd.read_csv(path, index_col=0)
    sub_df = get_condi_sub(df)
    ret = sub_df['real_return'].mean()
    rets[date] = ret

# %%
(pd.Series(rets).sort_index() + 1).cumprod().plot(figsize = (16, 8))
plt.show()

