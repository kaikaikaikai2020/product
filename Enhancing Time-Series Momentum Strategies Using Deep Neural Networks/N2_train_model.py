import torch
import pandas as pd
import os
import numpy as np
from torch.utils.data import DataLoader, Dataset
from tqdm import tqdm
from torch.utils.data.dataset import random_split, Subset
import matplotlib.pyplot as plt

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# device = torch.device('cpu')
torch.random.manual_seed(2020)

# %% read data
df = pd.read_csv('./feat_data/A.csv', index_col='tradeDate')
x_cols = ['rocp', 'orocp', 'lrocp', 'hrocp', 'ma', 'macp', 'dif', 'dea', 'macd',
          'difrocp', 'dearocp', 'macdrocp', 'rsi', 'rsirocp', 'vrocp', 'up',
          'mid', 'low', 'vmarocp', 'vmacp', 'rocpvrocp']
stk_codes = [i for i in os.listdir('./feat_data') if '.csv' in i]


# %% rolling standard scale
def roll_standard_scale(df_, x_cols=x_cols, window=50):
    df = df_.copy()
    for i in x_cols:
        mean_ = df[i].rolling(window).mean()
        std_ = df[i].rolling(window).std()

        df[i] = (df[i] - mean_) / std_

        # 做上下限处理
        def limit_val(x):
            if x > 2.6:
                x = 2.6
            if x < -2.6:
                x = -2.6
            return x

        df[i] = df[i].apply(limit_val)
    df['y'] = df['chgPct'] / df['sigma']
    df['eval_ret'] = df['nextPct']
    return df[x_cols + ['y', 'eval_ret']]


# %% leave 1500 data to test 1500 data to train
data_dict = {}


def read_data(stk_code):
    df = pd.read_csv('./feat_data/%s.csv' % stk_code, index_col='tradeDate')
    if df.dropna().shape[0] > 3000:  # 保证数据量
        df_scale = roll_standard_scale(df)
        data_dict[stk_code.split('.')[0]] = df_scale.dropna()


for stk_code in stk_codes: read_data(stk_code.split('.')[0])

# %% split train and test
arr_data = {}
for stk in data_dict.keys():
    df = data_dict[stk].copy()
    df_v = df.values
    df_test = df.loc['2017-*-*':, :]
    n, m, l = df.shape[0], df_test.shape[0], df.shape[1] - 2
    window = 63
    data_x = np.zeros(shape=(n - window, window, l))
    data_y = np.zeros(shape=(n - window, 1))

    for i in range(n - window):
        data_x[i] = df_v[i: i + window, :-2]
        data_y[i] = df_v[i + window, -2]
    print(stk, m, data_x[-m:].shape)
    arr_data[stk] = (data_x[:-m], data_y[:-m], data_x[-m:], df_test)


# %%
class tensor_data(Dataset):

    def __init__(self, x, y):
        self.x = x
        self.y = y.reshape(-1, 1)

    def __getitem__(self, item):
        x_i = torch.FloatTensor(self.x[item])
        y_i = torch.FloatTensor(self.y[item])
        return x_i, y_i

    def __len__(self):
        return len(self.x)


# %%
def train(data, model, loss_func, optimizer):
    train_loss = 0
    flag = 0
    data_loader = DataLoader(data, batch_size=256, shuffle=True)
    for (x_b, y_b) in tqdm(data_loader):
        optimizer.zero_grad()
        x_b = x_b.to(device)
        y_b = y_b.to(device)

        # model computer
        y_pred = model(x_b)
        loss = loss_func(y_pred, y_b)
        flag += 1

        train_loss += loss.item()
        loss.backward()
        optimizer.step()

    return train_loss / flag


def test(data, model, loss_func):
    train_loss = 0
    flag = 0
    batch_size = 256
    data_loader = DataLoader(data, batch_size=batch_size, shuffle=False)

    pred_arr = np.zeros(shape=(len(data), 1))

    for (x_b, y_b) in tqdm(data_loader):
        x_b = x_b.to(device)
        y_b = y_b.to(device)

        # model computer
        with torch.no_grad():
            y_pred = model(x_b)
            loss = loss_func(y_pred, y_b)
            train_loss += loss.item()
            # save pred
            if device.type == 'cuda':
                y_pred = y_pred.to('cpu')
            pred_arr[flag * batch_size: (flag + 1) * batch_size] = y_pred.numpy()

            flag += 1

    return train_loss / flag, pred_arr


# %%
train_x, train_y, test_x, eval_dfs = [], [], [], []
for k, v in arr_data.items():
    train_x.append(v[0])
    train_y.append(v[1])
    test_x.append(v[2])
    eval_dfs.append((k, v[3]))

train_data = tensor_data(np.vstack(train_x), np.vstack(train_y))
# 预测填充补齐，不参与计算
test_data = tensor_data(np.vstack(test_x), np.zeros(shape=(np.vstack(test_x).shape[0], 1)))
train_data, val_data = random_split(train_data, [int(len(train_data) * 0.7),
                                                 len(train_data) - int(len(train_data) * 0.7)])

# %%
from algs.lstm_model import *

loss_type = {'reg': reg_loss(), 'binary': binary_loss(), 'shp': sharpe_loss(), 'avg': avg_return_loss()}
for loss_type, loss_func in loss_type.items():
    print(loss_type)
    n_round = 400
    model = lstm_model(input_feat=21)
    model.to(device)
    loss_func.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    earlystop = EarlyStopping(patience=30, verbose=False, checkpoint='./model.pt')

    for i in range(n_round):
        loss = train(train_data, model, optimizer=optimizer, loss_func=loss_func)
        loss_v, _ = test(val_data, model, loss_func=loss_func)
        print(loss, loss_v)
        earlystop(loss_v, model)
        if earlystop.early_stop:
            # print('the model has done')
            break
    model.load_state_dict(torch.load('./model.pt'))
    _, pred = test(test_data, model, loss_func=loss_func)
    evals = []
    for i in eval_dfs:
        eval_df = i[1].copy()
        eval_df['code'] = i[0]
        evals.append(eval_df)
    eval_df = pd.concat(evals)
    eval_df['pred_val'] = pred
    if not os.path.exists('./eval_data'): os.mkdir('./eval_data')
    eval_df.to_csv('./eval_data/%s.csv' % loss_type)
