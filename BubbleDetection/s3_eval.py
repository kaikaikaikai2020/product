from glob import glob
import pandas as pd
import matplotlib.pyplot as plt

# %%
indu_list = glob('./signal/*.csv')

# %%
base_df = pd.read_csv(indu_list[1], index_col=0)

# %%
signals = []
for path in indu_list:
    base_df = pd.read_csv(path, index_col=0)
    base_df['chgpct'] = base_df['closeIndex'].shift(-1) / base_df['closeIndex'] - 1
    signal = (base_df['pred'] > 0.2)  # & ((base_df['openIndex'] - base_df['closeIndex'].shift(1)) > 0)
    s = base_df['chgpct']
    signals.append(s)

ret_df = pd.concat(signals, axis=1)
ret_s = ret_df.apply(lambda x: x[x != 0].mean(), axis=1).fillna(0)
(ret_s + 1).cumprod().plot(label = 'base')

# %%
signals = []
for path in indu_list:
    base_df = pd.read_csv(path, index_col=0)
    base_df['chgpct'] = base_df['closeIndex'].shift(-1) / base_df['closeIndex'] - 1
    signal = (base_df['pred'] > 0.2)  # & ((base_df['openIndex'] - base_df['closeIndex'].shift(1)) > 0)
    s = signal * base_df['chgpct']
    signals.append(s)

ret_df = pd.concat(signals, axis=1)
ret_s = ret_df.apply(lambda x: x[x != 0].mean(), axis=1).fillna(0)
(ret_s + 1).cumprod().plot(label = 'lambda 0.2')
# plt.legend()
# plt.show()

# %%
signals = []
for path in indu_list:
    base_df = pd.read_csv(path, index_col=0)
    base_df['chgpct'] = base_df['closeIndex'].shift(-1) / base_df['closeIndex'] - 1
    signal = (base_df['pred'] > 0.5) & ((base_df['openIndex'] - base_df['closeIndex'].shift(1)) > 0)
    s = signal * base_df['chgpct']
    signals.append(s)

ret_df = pd.concat(signals, axis=1)
ret_s = ret_df.apply(lambda x: x[x != 0].mean(), axis=1).fillna(0)
(ret_s + 1).cumprod().plot(label = 'lambda 0.5')
# plt.legend()
# plt.show()

# %%
signals = []
for path in indu_list:
    base_df = pd.read_csv(path, index_col=0)
    base_df['chgpct'] = base_df['closeIndex'].shift(-1) / base_df['closeIndex'] - 1
    signal = (base_df['pred'] > 0.8) & ((base_df['openIndex'] - base_df['closeIndex'].shift(1)) > 0)
    s = signal * base_df['chgpct']
    signals.append(s)

ret_df = pd.concat(signals, axis=1)
ret_s = ret_df.apply(lambda x: x[x != 0].mean(), axis=1).fillna(0)
(ret_s + 1).cumprod().plot(label = 'lambda 0.8')
plt.legend()
plt.show()
