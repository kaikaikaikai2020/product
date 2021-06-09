import pandas as pd
import numpy as np
import os
from tqdm import tqdm
from bubble_alg import bubble_boost

# %%
df = pd.read_csv('./data/indu_index.csv', index_col=0, parse_dates=['tradeDate'])
indu_id = df['industryID1'].unique().tolist()

# %%
for indu in indu_id:
    sub_df = df[df['industryID1'] == indu].copy()
    sub_df.index = sub_df.pop('tradeDate')
    sub_df.sort_index(inplace=True)

    # get weekday data
    sub_df['weekday'] = sub_df.index.map(lambda x: x.weekday())
    pre_df, base_df = sub_df[sub_df['weekday'] == 1], sub_df[sub_df['weekday'].shift(1) == 1]

    if base_df.shape[0] == pre_df.shape[0]:
        print('the shape is ok')

    else:
        print('the pre_df in more')
        pre_df = pre_df.iloc[:-1, ]
    
    # get price data
    base_c, pre_c = base_df['closeIndex'].values, pre_df['closeIndex'].values
    
    window = 150  # use 3 years data to get signal
    adf_quantile = np.zeros(shape=(base_c.shape[0] - window + 1,))
    flag = 0
    for i in tqdm(range(window, base_c.shape[0] + 1)):
        sub_pre = pre_c[i - window: i]
        model = bubble_boost(sub_pre)
        adfs = model.fast_fit()
        t_adf = model.get_true_adf()
        adf_quantile[flag] = (adfs < t_adf).mean()
        flag += 1
        # print((adfs < t_adf).mean(), t_adf)
    
    # saving signal
    base_df['pred'] = 0
    base_df['pred'].iloc[-len(adf_quantile):] = adf_quantile
    
    # saving data
    if not os.path.exists('./signal'): os.mkdir('./signal')
    base_df.to_csv('./signal/%s.csv' % indu)
    print('%s indu has done' % indu)
