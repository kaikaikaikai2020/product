import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt

# %%
df = pd.read_csv('./data/indu_week.csv', index_col=0)
indu_id = df['industryID1'].unique().tolist()
sub_df = df[df['industryID1'] == indu_id[0]].copy()
sub_df.index = sub_df.pop('endDate')
sub_df.sort_index(inplace=True)
sub_df = sub_df.iloc[4:]  # 对刚上的指数进行一些基本操作
y, delta_y = sub_df['closePrice'].values, sub_df['chg'].values

# %%
delta_y = y[1:] - y[:-1]
n = len(y)
delta_y, delta_ym1 = delta_y[1:], delta_y[:-1]
y, ym1, y1 = y[2:], y[1:-1], y[1]

# %%
x = sm.add_constant(np.vstack([ym1, delta_ym1]).T)
result = sm.OLS(delta_y, x).fit()
print(result.tvalues[1], 'stats model t-ratio')

# %%
np.linalg.inv(x.T @ x) @ x.T @ delta_y

# %%
xTx = np.linalg.inv(x.T @ x)
parmas = xTx @ x.T @ delta_y
e = delta_y - x @ parmas
print(parmas[1] / np.sqrt(xTx[1, 1] * e @ e / (1077 - 3)), 'numpy way')
#
# xTx[1, 1] *

# %%
n = x.shape[0]
x_tensor = np.zeros(shape=(10, x.shape[0], x.shape[1]))
y_tensor = np.zeros(shape=(10, y.shape[0], 1))
for i in range(10):
    x_tensor[i] = x
    y_tensor[i] = delta_y.reshape(-1, 1)

# %%
x_trans = x_tensor.transpose(0, 2, 1)
xTx = np.matmul(x_trans, x_tensor)
xTx_inv = np.linalg.inv(xTx)
params = np.matmul(np.matmul(xTx_inv, x_trans), y_tensor)
e_tensor = (y_tensor - np.matmul(x_tensor, params)).reshape(10, -1)
t_vals = params[:, 1, :].ravel() / np.sqrt((xTx_inv[:, 1, 1].ravel() * (e_tensor * e_tensor).sum(axis=1) / (n - 3)))
print(t_vals, 'tensor ways')

