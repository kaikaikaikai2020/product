import torch
from torch import nn as nn
import pandas as pd
from glob import glob
from tqdm import tqdm
from torch.utils.data import TensorDataset, DataLoader
import copy
import os

# %%
month_list = glob('./month_data/*.csv')
month_list.sort()

dfs = {}
for path in tqdm(month_list):
    date = path[-11:-4]
    df = pd.read_csv(path, index_col=0)
    df_ = df.iloc[:, 3:-3]
    dfs[date] = (df[['real_return']], df_)


# %%
def deep_model(input_size, drop_out=0.5):
    model = nn.Sequential(
        nn.Linear(input_size, 512),
        nn.Dropout(drop_out),
        nn.ReLU(),
        nn.Linear(512, 200),
        nn.Dropout(drop_out),
        nn.ReLU(),
        nn.Linear(200, 200),
        nn.Dropout(drop_out),
        nn.ReLU(),
        nn.Linear(200, 200),
        nn.Dropout(drop_out),
        nn.ReLU(),
        nn.Linear(200, 128),
        nn.Dropout(drop_out),
        nn.ReLU(),
        nn.Linear(128, 3),
        nn.Softmax(dim=1),
    )
    return model


def train_model(train_td):
    model.train()
    data_loader = DataLoader(train_td, shuffle=True, batch_size=256)
    losses, flag = 0, 0
    for x_i, y_i in tqdm(data_loader):
        y_pred = model(x_i)
        loss = loss_func(y_pred, y_i)
        opt.zero_grad()
        loss.backward()
        opt.step()
        losses += loss.item()
        flag += 1
    return losses / flag


def eval_model(test_x, test_y=None, model=None, is_test=False):
    model.eval()
    with torch.no_grad():
        pred_ = model(test_x)
        if not is_test:
            loss = loss_func(pred_, test_y)
            return loss.item()
        else:
            return (pred_[:, -1] - pred_[:, 0]).data.numpy()


# %%
month_step = 6

for i in range(len(month_list) - month_step):
    train_dates = month_list[i: i + month_step - 1]
    train_data = []
    for path in train_dates:
        date = path[-11:-4]
        _, df = dfs[date]
        train_data.append((torch.tensor(df.iloc[:, :-1].values, dtype=torch.float32),
                           torch.tensor(pd.qcut(df['return'], 3, range(3)).values,
                                        dtype=torch.long)))

    train_td = TensorDataset(torch.cat([i[0] for i in train_data], axis=0),
                             torch.cat([i[1] for i in train_data], axis=0))

    val_date = month_list[i + month_step - 1][-11:-4]
    val_df, val_data = dfs[val_date]
    val_x = torch.tensor(val_data.iloc[:, :-1].values, dtype=torch.float32)
    val_y = torch.tensor(pd.qcut(val_data['return'], 3, range(3)).values,
                         dtype=torch.long)

    test_date = month_list[i + month_step][-11:-4]
    test_df, test_data = dfs[test_date]
    test_x = torch.tensor(test_data.iloc[:, :-1].values, dtype=torch.float32)

    model = deep_model(244)
    loss_func = nn.CrossEntropyLoss()
    opt = torch.optim.Adam(model.parameters(), lr=1e-4)

    p_flag = 0
    patience = 10
    val_best = 1000
    for i in range(100):
        train_loss = train_model(train_td)
        val_loss = eval_model(val_x, val_y, model)
        print('train loss: %.4f | val loss: %.4f ' % (train_loss, val_loss))
        if val_loss < val_best:
            val_best = val_loss
            best_model = copy.deepcopy(model)
            p_flag = 0
        p_flag += 1
        if p_flag > 10:
            break
    pred_ = eval_model(test_x, model=model, is_test=True)
    test_df['pred'] = pred_
    if not os.path.exists('./deep_result/'): os.mkdir('./deep_result/')
    test_df.to_csv('./deep_result/%s.csv' % test_date)
    print('========== %s has done' % test_date)
