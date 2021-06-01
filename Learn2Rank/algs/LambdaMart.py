import xgboost as xgb
from xgboost import DMatrix
import numpy as np

# %%
ws = np.random.randn(16, 1)
x_train_1 = np.random.randn(1000, 16)
x_train_2 = np.random.randn(1000, 16)
x_val = np.random.randn(1000, 16)
y_train = np.exp((np.vstack([x_train_1, x_train_2]) @ ws).ravel())
# y_train[y_train.argsort()] = np.arange(y_train.shape[0])
y_val = np.exp((x_val @ ws).ravel())
# y_val[y_val.argsort()] = np.arange(y_val.shape[0])
train_M = DMatrix(np.vstack([x_train_1, x_train_2]), y_train)
val_M = DMatrix(x_val, y_val)
group_x = np.array([len(x_train_1) * 2])
train_M.set_group(group_x)
# %%
# params = {'objective': 'rank:pairwise', 'eta': 0.01, 'eval_metric': 'ndcg@1000'}
# model = xgb.train(params, train_M, evals=[(train_M, 'train'), (val_M, 'val')], num_boost_round=200)

# %%
model = xgb.XGBRanker(
    booster='gbtree',
    objective='rank:pairwise',
    eval_metric='ndcg@100',
    random_state=42,
    learning_rate=0.1,
    colsample_bytree=0.9,
    eta=0.05,
    max_depth=6,
    n_estimators=110,
    subsample=0.75
)

model.fit(np.vstack([x_train_1, x_train_2]), y_train, eval_set=[(x_val, y_val)],
          eval_group=[np.array([len(y_val)])],
          group=group_x, verbose=True, early_stopping_rounds= 10)
