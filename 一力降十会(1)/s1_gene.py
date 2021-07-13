import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from algs.genetic import delay_rolling_func, cross_delay_func, single_func, cross_func, \
    generate_data_by_sign, generate_sign, generate_norepeat_signs, correlation_filter
np.random.seed(2021)

# %%
df = pd.read_csv('./data/indu_index.csv', index_col=0)
indu_code = df['industryID1'].unique().tolist()

indu_dfs = []
for code in indu_code:
    sub_df = df[df['industryID1'] == code].sort_values('tradeDate')
    indu_dfs.append(sub_df)


# %%
def get_base_data(sub_df):
    sub_df['dif_close'] = sub_df['closeIndex'].diff() / sub_df['closeIndex'].shift(1)
    sub_df['dif_open'] = sub_df['openIndex'].diff() / sub_df['openIndex'].shift(1)
    sub_df['dif_high'] = sub_df['highestIndex'].diff() / sub_df['highestIndex'].shift(1)
    sub_df['dif_low'] = sub_df['lowestIndex'].diff() / sub_df['lowestIndex'].shift(1)
    sub_df['dif_vol'] = sub_df['turnoverVol'].diff() / sub_df['turnoverVol'].shift(1)

    df_x = sub_df[['closeIndex', 'openIndex', 'highestIndex', 'lowestIndex',
                   'turnoverVol', 'turnoverValue',
                   'dif_close', 'dif_open', 'dif_high', 'dif_low', 'dif_vol', 'CHG']]
    df_y = sub_df['CHGPct'].shift(-1)
    return (df_x - df_x.mean()) / df_x.std(), df_y

def stable_filter(df_x, df_y, chs_number=100, split_part=10):
    corrs = []
    df_xcor = df_x.dropna()
    df_ycor = df_y.loc[df_xcor.index]
    n = df_x.shape[0] // split_part + 1
    for i in range(split_part):
        sub_x, sub_y = df_xcor.iloc[i * n: (i + 1) * n], df_ycor.iloc[i * n: (i + 1) * n]
        corr_df = sub_x.corrwith(sub_y, method='spearman')
        corrs.append(corr_df)
    corr_df = pd.concat(corrs, axis=1)

    # ic = corr_df.mean(axis=1) / corr_df.std(axis=1)
    # ic = (corr_df * np.linspace(0, 1, split_part)).sum(axis = 1)
    ic = corr_df.sum(axis=1)
    # ic = corr_df.sum(axis=1) * corr_df.applymap(np.sign).mean(axis=1)

    sort_cols = ic.abs().sort_values(ascending=False)
    df_sort = df_x[sort_cols.index]
    corr_df = df_sort.corr()
    corr_df.replace(1, 0, inplace=True)
    cols_continue = sort_cols.index.tolist()
    flag = 0
    for col in cols_continue:
        flag += 1
        sub = corr_df[corr_df[col].abs() > 0.6]
        if len(sub.index) != 0:
            for sub_col in sub.index:
                cols_continue.remove(sub_col)
        corr_df = corr_df.loc[cols_continue, cols_continue]
        if flag > chs_number:
            break
    return corr_df.index.tolist()[:chs_number]

# %%
test_size = 1000

df_x, df_y = get_base_data(sub_df)
for i in range(4):
    sign_list = generate_norepeat_signs(df_x.columns, 300)
    sub_xs = [generate_data_by_sign(df=df_x, signs=signs) for signs in sign_list]
    sub_dfx = pd.concat(sub_xs, axis=1)
    sub_dfx = sub_dfx.loc[:, ~sub_dfx.columns.duplicated()]

    df_x = pd.concat([df_x, sub_dfx], axis=1)
    df_x = df_x.loc[:, ~df_x.columns.duplicated()]
    # avoid nan value
    df_x = df_x[df_x.columns[df_x.isna().mean() < 0.2]]
    # df_xtrain, df_xtest, df_ytrain, df_ytest = df_x.iloc[:-test_size], df_x.iloc[-test_size:], \
    #                                            df_y.iloc[:-test_size], df_y.iloc[-test_size:]

    # chs_cols = correlation_filter(df_x, df_y)
    chs_cols = stable_filter(df_x, df_y)

    # ic = df_x.corrwith(df_y).drop_duplicates().dropna()
    # # get ic abs top
    # chs_cols = ic.abs().sort_values(ascending=False).iloc[:100]
    df_x = df_x[chs_cols]
    cor = df_x.corrwith(df_y, method='spearman').drop_duplicates().dropna()
    print(cor)

# %%
sign = np.zeros_like(df_y)
if cor.iloc[0] < 0:
    sign_bool = df_x.iloc[:, 0] < df_x.iloc[:, 0].quantile(0.5)
else:
    sign_bool = df_x.iloc[:, 0] > df_x.iloc[:, 0].quantile(0.5)
sign[sign_bool] = 1
(df_y + 1).cumprod().plot()
(df_y * sign + 1).cumprod().plot()
plt.show()


# %%
#
(df_y.iloc[-test_size:] + 1).cumprod().plot()
((df_y * sign).iloc[-test_size:] + 1).cumprod().plot()
plt.show()
