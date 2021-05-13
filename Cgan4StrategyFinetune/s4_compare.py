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

x_train, x_test = x_data[:c], x_data[c:]
y_train, y_test = y_data[:c], y_data[c:]
y_train[y_train > 0] = 1
y_train[y_train <= 0] = 0
# %%
test_fold = np.full((y_train.shape[0]), 0, dtype=int)
test_fold[:3000] = -1

# %% gbdt case
from sklearn.ensemble import GradientBoostingRegressor, GradientBoostingClassifier
from sklearn.model_selection import GridSearchCV, PredefinedSplit
import joblib

ps = PredefinedSplit(test_fold)

model = GradientBoostingRegressor(random_state=2021)
params = {'n_estimators': [20, 50, 100, 200],
          'learning_rate': [0.001, 0.01, 0.1, 1],
          'max_depth': [1, 3, 5]}
gsearch = GridSearchCV(estimator=GradientBoostingClassifier(random_state=2021),
                       param_grid=params, cv=ps,
                       scoring='accuracy', n_jobs=-1)
gsearch.fit(x_train, y_train.ravel())

# %%
cgan_model = joblib.load('./cgan_model.pkl').best_estimator_
base_model = gsearch.best_estimator_

# %%
plt.figure(figsize=(16, 8))
plt.plot((y_test.ravel() + 1).cumprod(), label='index')
plt.plot((cgan_model.predict(x_test) * y_test.ravel() + 1).cumprod(), label='cgan')
plt.plot((base_model.predict(x_test) * y_test.ravel() + 1).cumprod(), label='base')
plt.legend()
plt.show()
