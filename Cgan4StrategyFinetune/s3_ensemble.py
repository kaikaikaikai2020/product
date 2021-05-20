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


# %% gbdt case
from sklearn.ensemble import GradientBoostingRegressor, GradientBoostingClassifier
from sklearn.model_selection import GridSearchCV, PredefinedSplit
import joblib
from tqdm import tqdm

best_models = []
for i in tqdm(range(50)):
    x, y = cgan_sim()
    y[y > 0] = 1
    y[y <= 0] = 0
    # split train test
    test_fold = np.full((y.shape[0]), 0, dtype=int)
    test_fold[:3000] = -1
    ps = PredefinedSplit(test_fold)

    model = GradientBoostingRegressor(random_state=2021)
    params = {'n_estimators': [20, 50, 100, 200],
              'learning_rate': [0.001, 0.01, 0.1, 1],
              'max_depth': [1, 3, 5]}
    gsearch = GridSearchCV(estimator=GradientBoostingClassifier(random_state=2021),
                           param_grid=params, cv=ps,
                           scoring='accuracy', n_jobs=-1)
    gsearch.fit(x, y.ravel())
    best_models.append(gsearch)
# %%
joblib.dump(best_models, 'cgan_ensemble_model.pkl')
