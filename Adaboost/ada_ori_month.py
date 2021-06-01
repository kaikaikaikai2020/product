import pandas as pd
import numpy as np
from algs.adaboost import *
import os
from scipy.stats import spearmanr
from tqdm import tqdm
import matplotlib.pyplot as plt
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier

data_list = os.listdir('./month_data')
data_list = [i for i in data_list if 'csv' in i]
data_list.sort()


# %%
def get_train_test_path(begin_mon='2010-01', window=12, base_dir='./month_data'):
    for begin_path in data_list:
        if begin_mon in begin_path: break
    beg_ind = data_list.index(begin_path)
    train_ind = data_list[beg_ind: beg_ind + window]
    test_ind = data_list[beg_ind + window]
    return [os.path.join(base_dir, i) for i in train_ind], os.path.join(base_dir, test_ind)


def processing_data_order(s):
    s_order = s.argsort()
    return s_order / s_order.max()


def processing_data(df, is_train=True):
    order_cols = df.columns[3: -3]
    df_order = df[order_cols]
    df_x, df_y = df_order.iloc[:, :-1], df_order.iloc[:, -1]
    if is_train:
        chs_up = df_y > np.quantile(df_y, 0.5)
        chs_down = df_y < np.quantile(df_y, 0.5)
        df_y[chs_up] = 1
        df_y[chs_down] = 0
        df_y = df_y[chs_up | chs_down]
        df_x = df_x[chs_up | chs_down].apply(lambda x: processing_data_order(x))
        return df_x, pd.DataFrame(df_y)
    else:
        df_x = df_x.apply(lambda x: processing_data_order(x))
        return df_x, df.loc[:, ['ticker', 'tradeDate', 'real_return']]


def get_train_test_data(train_paths, test_path):
    train_datas = []
    for train_path in train_paths:
        df_train = pd.read_csv(train_path, index_col=0)
        df_train_x, df_train_y = processing_data(df_train, is_train=True)
        train_datas.append((df_train_x, df_train_y))
    df_test = pd.read_csv(test_path, index_col=0)
    df_test.dropna(inplace = True)
    df_test, df_ret = processing_data(df_test, is_train=False)
    x_train = pd.concat([i[0] for i in train_datas])
    y_train = pd.concat([i[1] for i in train_datas])
    return (x_train, y_train), df_test, df_ret



# %%
window = 12
ic_dict = {}
res_dict = {}
for date in tqdm(data_list[:-window]):
    train_paths, test_path = get_train_test_path(begin_mon= date, window= window)
    test_date = test_path.split('/')[-1].split('.')[0]
    (x_train, y_train), x_test, df_ret = get_train_test_data(train_paths, test_path)
    # model = adaboost_weak_cls(n_models= 10, Q= 5)
    model = AdaBoostClassifier(DecisionTreeClassifier(max_depth=2), n_estimators=30, learning_rate=0.8)
    model.fit(x_train.values, y_train.values.ravel())
    train_pred = model.predict_proba(x_train)[:, 1]
    test_pred = model.predict_proba(x_test)[:, 1]
    train_ic = spearmanr(train_pred.ravel(), y_train.values.ravel())[0]
    test_ic = spearmanr(test_pred.ravel(), df_ret['real_return'].values)[0]
    df_ret['els'] = test_pred
    ic_dict[test_date] = test_ic
    res_dict[test_date] = df_ret
    # print(train_ic, test_ic)


# %%
ret_dfs = []
n = 5
for k in res_dict.keys():
    try:
        res_df = res_dict[k]
        ret_df = res_df.groupby(pd.qcut(res_df['els'], n, range(n)))['real_return'].mean()
        ret_df.name = k.split('\\')[1]
        ret_dfs.append(ret_df)
    except :
        continue


ret_df = pd.concat(ret_dfs, axis = 1).T

# %%
(ret_df + 1).cumprod().plot()
(ret_df[n - 1] - ret_df[0] + 1).cumprod().plot()
plt.title('adaboost_ori')
plt.show()