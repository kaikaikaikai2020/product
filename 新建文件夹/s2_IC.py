from glob import glob
import pandas as pd
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt

# %%
month_list = glob('./month_data/*.csv')
month_list.sort()

dfs = {}
for path in tqdm(month_list):
    date = path[-11:-4]
    df = pd.read_csv(path, index_col=0)
    df_ = df.iloc[:, 3:-3]
    df['month'] = date
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
    if sig > 0:
        test_df.sort_values(fac_id, ascending=False, inplace=True)
    else:
        test_df.sort_values(fac_id, ascending=True, inplace=True)
    ret = test_df['real_return'].iloc[:35].mean() # - test_df['real_return'].iloc[-35:].mean()
    rets[test_date] = ret

# %% max ic
(pd.Series(rets).sort_index() + 1 - 0.005).cumprod().plot(figsize=(16, 8))
plt.show()
