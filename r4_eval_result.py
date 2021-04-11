import pandas as pd
import matplotlib.pyplot as plt

# %%
df = pd.read_csv('./moment.csv', index_col=0)
df.index = df.index.astype(str)

# %%
month_list = df.index.unique().tolist()
month_list.sort()


# %%
def read_index():
    df = pd.read_csv('./index_data/sh_index.csv', index_col=0)

    def split_date(x):
        y, m, d = x.split('-')
        return y + m + d

    df['trade_date'] = df['tradeDate'].apply(split_date)
    df['month'] = df['trade_date'].apply(lambda x: x[:6])
    return df.groupby('month')['CHGPct'].apply(lambda x: (x + 1).prod() - 1)


# %%
n = 10
sub_ms = []
sub_rs = []
for month in month_list:
    sub_df = df.loc[month]
    sub_r = sub_df.groupby(pd.qcut(sub_df['mom_residual'], n, range(1, n + 1)))['chg_pct'].mean()
    sub_m = sub_df.groupby(pd.qcut(sub_df['mom_convent'], n, range(1, n + 1)))['chg_pct'].mean()
    sub_r.name, sub_m.name = month, month
    sub_ms.append(sub_m)
    sub_rs.append(sub_r)

# %%
residual_df = pd.concat(sub_rs, axis=1).T
convent_df = pd.concat(sub_ms, axis=1).T
df_index = read_index()
# %%
(df_index.loc[residual_df.index] + 1).cumprod().plot(figsize=(16, 8), label='index')
(residual_df[n] - residual_df[1] + 1).cumprod().plot(label='Residual')
(convent_df[n] - convent_df[1] + 1).cumprod().plot(label='conventional')
plt.legend()
plt.show()
