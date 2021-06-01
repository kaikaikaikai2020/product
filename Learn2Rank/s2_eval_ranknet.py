import pandas as pd
import numpy as np
from glob import glob
import matplotlib.pyplot as plt

# %%
res_list = glob('./ranknet/*.csv')
res_list.sort()


# %%
def read_stk(path):
    df = pd.read_csv(path)
    # df = df[df['next_month_ret'].abs() < 1]
    return df


# %%
ret_dict = {}
for path in res_list:
    df = read_stk(path)
    date = path[-11:-4]
    if date not in ret_dict:
        ret_dict[date] = {}
    for col in df.columns:
        if 'pred' in col:
            sort_ind = df[col].sort_values(ascending=False).index
            top_bot = df.loc[sort_ind[:100], 'real_return'].apply(lambda x: -1 if x < -1 else x).mean() - \
                      df.loc[sort_ind[-100:], 'real_return'].apply(lambda x: 1 if x > 1 else x).mean()
            # top_bot = df.loc[sort_ind[:100], 'next_month_ret'].apply(lambda x: -1 if x < -1 else x).mean()
            ret_dict[date][col] = top_bot

# %%
ret_df = pd.DataFrame(ret_dict).T.sort_index()
(ret_df + 1).cumprod().plot(figsize=(16, 8))
plt.show()

# %%
print(ret_df.mean() / ret_df.std() * np.sqrt(12))
