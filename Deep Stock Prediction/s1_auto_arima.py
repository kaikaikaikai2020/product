import pandas as pd
import numpy as np
import pmdarima as pm
import os
import matplotlib.pyplot as plt
from tqdm import tqdm

# %%
data_dict = {}
file_list = os.listdir('./stk_data')
for file_name in file_list:
    if 'csv' in file_name and 'b' not in file_name:
        print(file_name)
        df = pd.read_csv('./stk_data/%s' % file_name, index_col=0)
        df.dropna(inplace=True)
        df['mid'] = (df['HIGH'] + df['LOW']) / 2
        data_dict[file_name.split('.')[0]] = df


# %%
def pred_step_by_step(test_data, model):
    test_size = len(test_data)
    preds = np.zeros(test_size)
    for i in tqdm(range(test_size)):
        new_pred = model.predict(n_periods=1)
        model.update([test_data[i]])
        preds[i] = new_pred
    return preds


# %%
for k, df in data_dict.items():
    train_data, test_data, test_df = df['mid'].iloc[:-1000], df['mid'].iloc[-1000:], df.iloc[-1000:].copy()

    # auto-arima model
    model = pm.auto_arima(train_data.values, stationary=False, trace=True, seasonal=True, m=12)
    # forecasts = model.predict(test.shape[0])
    # plt.plot(model.predict_in_sample()[1:])
    # plt.plot(train_data.values[1:])
    # plt.show()

    # pred step by step
    preds = pred_step_by_step(test_data.values, model)
    test_df['pred'] = preds
    test_df['mid'].plot()
    test_df['pred'].plot()
    plt.show()
    if not os.path.exists('./pred_data'): os.mkdir('./pred_data')
    test_df.to_csv('./pred_data/%s.csv' % k)
