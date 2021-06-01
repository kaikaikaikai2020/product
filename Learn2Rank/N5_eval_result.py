import pandas as pd
import numpy as np
from glob import glob
import matplotlib.pyplot as plt


# %%
def read_stk(path):
    df = pd.read_csv(path)
    return df


def get_result(paths_dir='./ranknet/*.csv'):
    res_list = glob(paths_dir)
    res_list.sort()
    ret_dict = {}
    for path in res_list:
        df = read_stk(path)
        date = path[-11:-4]
        if date not in ret_dict:
            ret_dict[date] = {}
        sort_ind = df['pred'].sort_values(ascending=False).index
        top_bot = df.loc[sort_ind[:100], 'real_return'].apply(lambda x: -1 if x < -1 else x).mean() - \
                  df.loc[sort_ind[-100:], 'real_return'].apply(lambda x: 1 if x > 1 else x).mean()

        ret_dict[date] = top_bot
    ret_df = pd.Series(ret_dict).sort_index()
    return ret_df


# %%
ret_rank = get_result('./ranknet/*.csv')
ret_list = get_result('./listmle/*.csv')
ret_lamb = get_result('./lambda/*.csv')

# %%
(ret_rank + 1).cumprod().plot(label='ranknet', figsize=(16, 8))
(ret_list + 1).cumprod().plot(label='listmle', figsize=(16, 8))
(ret_lamb + 1).cumprod().plot(label='lambdamart', figsize=(16, 8))
plt.legend()
plt.show()

# %%
def get_adj_result(paths_dir='./ranknet/*.csv'):
    res_list = glob(paths_dir)
    res_list.sort()
    ret_dict = {}
    for path in res_list:
        df = read_stk(path)
        date = path[-11:-4]
        if date not in ret_dict:
            ret_dict[date] = {}
        sort_ind = df['pred'].sort_values(ascending=False).index
        top_df = df.loc[sort_ind[:100], ['real_return', 'market_value']]
        top = (top_df['real_return'] * top_df['market_value'].apply(np.sqrt) / top_df['market_value'].apply(
            np.sqrt).sum()).sum()
        bot_df = df.loc[sort_ind[-100:], ['real_return', 'market_value']]
        bot = (bot_df['real_return'] * bot_df['market_value'].apply(np.sqrt) / bot_df['market_value'].apply(
            np.sqrt).sum()).sum()

        ret_dict[date] = top - bot
    ret_df = pd.Series(ret_dict).sort_index()
    return ret_df


# %%
ret_rank = get_adj_result('./ranknet/*.csv')
ret_list = get_adj_result('./listmle/*.csv')
ret_lamb = get_adj_result('./lambda/*.csv')

# %%
(ret_rank + 1).cumprod().plot(label='ranknet', figsize=(16, 8))
(ret_list + 1).cumprod().plot(label='listmle', figsize=(16, 8))
(ret_lamb + 1).cumprod().plot(label='lambdamart', figsize=(16, 8))
plt.legend()
plt.show()
