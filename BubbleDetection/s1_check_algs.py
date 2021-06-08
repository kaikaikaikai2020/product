import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# %%
df = pd.read_csv('./data/indu_index.csv', index_col=0, parse_dates=['tradeDate'])
indu_id = df['industryID1'].unique().tolist()

# %%
sub_df = df[df['industryID1'] == indu_id[0]].copy()
sub_df.index = sub_df.pop('tradeDate')
sub_df.sort_index(inplace=True)

# %%
sub_df['weekday'] = sub_df.index.map(lambda x: x.weekday())
base_df, pre_df = sub_df[sub_df['weekday'] == 1], sub_df[sub_df['weekday'].shift(1) == 1]


# %%
from bubble_alg import bubble_boost

base_c, pre_c = base_df['closeIndex'].values, pre_df['closeIndex'].values
window = 150  # 大约3年
adf_quantile = np.zeros(shape=(base_c.shape[0] - window + 1,))
flag = 0
for i in tqdm(range(window, base_c.shape[0] + 1)):
    sub_pre = pre_c[i - window: i]
    model = bubble_boost(sub_pre)
    adfs = model.fast_fit()
    t_adf = model.get_true_adf()
    adf_quantile[flag] = (adfs < t_adf).mean()
    flag += 1
    print((adfs < t_adf).mean(), t_adf)

# %%
base_df['pred'] = 0
base_df['pred'].iloc[-len(adf_quantile):] = adf_quantile

# %%
base_df['closeIndex'].plot()
base_df[base_df['pred'] > 0.6]['closeIndex'].plot()
plt.show()

# %%
base_df['chgpct'] = base_df['closeIndex'].shift(-1) / base_df['closeIndex'] - 1
signal = (base_df['pred'] > 0.5)  # & (base_df['CHGPct'] > 0)
(signal * base_df['chgpct'] + 1).cumprod().plot(label = 'lambda = 0.5')
signal = (base_df['pred'] > 0.8)  # & (base_df['CHGPct'] > 0)
(signal * base_df['chgpct'] + 1).cumprod().plot(label = 'lambda = 0.8')
signal = (base_df['pred'] > 0.2)  # & (base_df['CHGPct'] > 0)
(signal * base_df['chgpct'] + 1).cumprod().plot(label = 'lambda = 0.2')
plt.legend()
plt.show()

# %%
