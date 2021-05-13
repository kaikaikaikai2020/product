import pandas as pd
import torch
import numpy as np
import matplotlib.pyplot as plt

# %%
g_model = torch.load('./g_model.pkl')
data = pd.read_csv('./data/hs_300.csv')

# %%
s = data['CHGPct'].values
window_size, n = 252, data.shape[0]
c = int((n - window_size) * .8)
x_data = np.array([s[i: i + window_size] for i in range(n - window_size)])
y_data = np.array([s[i + window_size] for i in range(n - window_size)]).reshape((-1, 1))


# %% simulate
def cgan_sim():
    x = torch.tensor(x_data[:c], dtype=torch.float32)
    z = torch.randn(size=(c, 252))
    return x.numpy(), g_model(torch.cat([z, x], axis=1)).data.numpy()


# %% generate training and val data
x_np_data, y_np_data = [], []
for i in range(1):
    x, y = cgan_sim()
    x_np_data.append(x)
    y_np_data.append(y)

x_np_data, y_np_data = np.vstack(x_np_data), np.vstack(y_np_data)

# %%
test_fold = np.full((y_np_data.shape[0]), 0, dtype=int)
for i in range(1):
    test_fold[i*x.shape[0]: i*x.shape[0] + 3000] = -1

# %% gbdt case
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import GridSearchCV, PredefinedSplit
import joblib

ps = PredefinedSplit(test_fold)

model = GradientBoostingRegressor(random_state=2021)
params = {'n_estimators': [20, 50, 100, 200],
          'learning_rate': [0.001, 0.01, 0.1, 1],
          'max_depth': [1, 3, 5]}
gsearch = GridSearchCV(estimator=GradientBoostingRegressor(random_state=2021),
                       param_grid=params, cv=ps,
                       scoring='neg_mean_squared_error', n_jobs=-1)


# %%
gsearch.fit(x_np_data, y_np_data.ravel())
joblib.dump(gsearch, 'cgan_reg_model.pkl')
# %%
