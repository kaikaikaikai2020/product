import pandas as pd
import numpy as np
import copy
import statsmodels.api as sm
from sklearn import linear_model
from sklearn import ensemble
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
from algs.adaboost import *
from scipy.stats import spearmanr

# %%
data_list = os.listdir('./month_data')
data_list = [i for i in data_list if 'csv' in i]
data_list.sort()


def get_train_test_path(begin_mon='2010-01', window=12, base_dir='./month_data'):
    for begin_path in data_list:
        if begin_mon in begin_path: break
    beg_ind = data_list.index(begin_path)
    train_ind = data_list[beg_ind: beg_ind + window]
    test_ind = data_list[beg_ind + window]
    return [os.path.join(base_dir, i) for i in train_ind], os.path.join(base_dir, test_ind)


def read_fac_data(begin_mon='2010-01', window=12, base_dir='./month_data'):
    train_paths, test_path = get_train_test_path(begin_mon, window, base_dir)
    # test_date = test_path.split('/')[-1].split('.')[0] # mac
    test_date = test_path.split('\\')[-1].split('.')[0]  # windows
    train_data = [pd.read_csv(path, index_col=0) for path in train_paths]
    test_data = pd.read_csv(test_path, index_col=0)
    return train_data, test_data, test_date


def processing_data_order(s):
    s_order = s.argsort()
    return s_order / s_order.max()


def processing_data(df, order_cols, y_col):
    df_x, df_y = df.loc[:, order_cols].copy(), df.loc[:, y_col].copy()
    chs_up = df_y > np.quantile(df_y, 0.7)
    chs_down = df_y < np.quantile(df_y, 0.3)
    df_y[chs_up] = 1
    df_y[chs_down] = -1
    df_y = df_y[chs_up | chs_down]
    df_x = df_x.apply(lambda x: processing_data_order(x))
    return pd.concat([df_x[chs_up | chs_down], pd.DataFrame(df_y)], axis=1)


def filter_icir(train_data, fac_cols, y_col):
    ic_df = pd.concat([data[fac_cols].corrwith(data[y_col]) for data in train_data], axis=1)
    icir = ic_df.mean(axis=1) / ic_df.std(axis=1)
    filter_cols = icir.abs().sort_values(ascending=False).iloc[:70].index
    return filter_cols.tolist()


def filter_by_adb(train_data, fac_cols, y_col):
    adb_df = pd.concat([processing_data(data, fac_cols, y_col) for data in train_data])
    adb_x, adb_y = adb_df[fac_cols].values, adb_df[[y_col]].values
    model = adaboost_weak_cls(n_models=30, Q=5)
    model.fit(adb_x, adb_y)
    filter_cols = [fac_cols[i[0]] for i in model.pred_list]
    return filter_cols


def main_gbdt_model(train_data, test_data, filter_cols, y_col):
    model = ensemble.GradientBoostingClassifier(n_estimators= 10)
    fil_df = pd.concat([data.loc[:, filter_cols + ['return']] for data in train_data], axis=0)
    result = model.fit(fil_df[filter_cols], fil_df[y_col].apply(lambda x: 1 if x > 0 else 0))
    pred = result.predict_proba(test_data[filter_cols])[:, 1]
    test_data['pred'] = pred
    return test_data[['real_return', 'return', 'pred']]


# %%
window = 12
# res_dict = {}
for mon_str in tqdm(data_list[:-window]):
    train_data, test_data, test_date = read_fac_data(begin_mon=mon_str, window=window)
    fac_cols, y_col = test_data.columns[3:-4].tolist(), 'return'
    filter_icir_cols = filter_icir(train_data, fac_cols, y_col)  # 利用 iric 进行因子筛选
    filter_adb = filter_by_adb(train_data, filter_icir_cols, y_col)  # 利用 adb 进行因子筛选
    test_result = main_gbdt_model(train_data, test_data, list(set(filter_adb)), y_col)  # 使用 adaboost 对因子整合
    # res_dict[test_date] = test_result # 保存结果
    if not os.path.exists('./adb_data'): os.mkdir('./adb_data')
    test_result.to_csv('./adb_data/res_%s.csv' % test_date)
    print(test_date, 'has done')
del train_data
del test_data

# %%
n = 5
ret_dfs = []
for mon_str in data_list[window:]:
    res_df = pd.read_csv('./adb_data/res_%s' % mon_str, index_col=0)
    ret_ser = res_df.groupby(pd.qcut(res_df['pred'], n, range(0, 5)))['real_return'].mean()
    ret_ser.name = mon_str
    ret_dfs.append(ret_ser)

ret_df = pd.concat(ret_dfs, axis=1).T
ret_df.sort_index(inplace=True)

# 因子分层检验，分层明显，说明因子是有效的
(ret_df + 1).cumprod().plot()
(ret_df[n-1] - ret_df[0] + 1).cumprod().plot()
plt.title('gbdt_model')
plt.show()

# %%
# 每月做多 分数最高的30只股票，做空 分数最低的30只股票
n_stocks = 30
rets = {}
for mon_str in data_list[window:]:
    res_df = pd.read_csv('./adb_data/res_%s' % mon_str, index_col=0)
    ret_short = res_df.sort_values('pred', ascending= True)['real_return'].iloc[:n].mean()
    ret_long = res_df.sort_values('pred', ascending= False)['real_return'].iloc[:n].mean()
    rets[mon_str.split('.')[0]] = ret_long - ret_short

(pd.Series(rets) + 1 - 0.005).cumprod().plot()
plt.show()
