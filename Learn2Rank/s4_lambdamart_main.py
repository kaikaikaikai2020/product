import xgboost as xgb
from xgboost import DMatrix
import pandas as pd
import numpy as np
from glob import glob
import os

# %%
month_list = glob('./month/*.csv')
month_list.sort()


def read_month_data(path):
    df = pd.read_csv(path, dtype={'ts_code': str}, index_col='ts_code')
    return (df.values[:, :-1], df.values[:, -1], df)


def get_train_val_test_data(month_path, train_window=3):
    begin_ind = month_list.index(month_path)
    train_paths = month_list[begin_ind: begin_ind + train_window]
    val_path = month_list[begin_ind + train_window]
    test_path = month_list[begin_ind + train_window + 1]
    return [read_month_data(p) for p in train_paths], \
           read_month_data(val_path), \
           read_month_data(test_path)

def quantile_x(x):
    x = x.ravel().copy()
    return np.exp(x)  # no neg score

# %%
train_window = 6
if not os.path.exists('./lambda_ed2'): os.mkdir('./lambda_ed2')
for ind in range(len(month_list) - train_window - 1):
    date = test_path = month_list[ind + train_window][-10:]
    train_data, val_data, test_data = get_train_val_test_data(month_list[ind], train_window=train_window)

    train_x = np.vstack([i[0] for i in train_data])
    train_y = np.hstack([quantile_x(i[1]) for i in train_data])

    val_x, val_y = val_data[0], quantile_x(val_data[1])
    model = xgb.XGBRanker(
        booster='gbtree',
        objective='rank:pairwise',
        eval_metric='ndcg@100',
        random_state=42,
        learning_rate=0.001,
        colsample_bytree=0.9,
        eta=0.05,
        max_depth=6,
        n_estimators=110,
        subsample=0.75
    )

    model.fit(train_x, train_y, eval_set=[(val_x, val_y)],
              eval_group=[np.array([len(val_y)])],
              group=np.array([len(train_y)]), verbose=False, early_stopping_rounds=10)

    pred_ = model.predict(test_data[0], ntree_limit=model.best_iteration)
    inds = pred_.argsort()
    test_data[2]['pred'] = pred_
    test_data[2].to_csv('./lambda_ed2/%s' % date)
    print('========== %s has done the ret is %.4f'% (date,
                                                     (test_data[1][inds[-100:]] - test_data[1][inds[:100]]).mean()))
