import pandas as pd
import numpy as np
import torch
from torch import nn
from torch.utils.data import DataLoader, random_split
from algs.ranknet import RankNet, pairwise_data
from glob import glob
from tqdm import tqdm
import os

torch.random.manual_seed(2021)

device = torch.device

# %%
month_list = glob('./month/*.csv')
month_list.sort()


# %%
def read_month_data(path, is_train=False):
    df = pd.read_csv(path, dtype={'ts_code': str}, index_col='ts_code')
    if is_train:
        return torch.tensor(df.values[:, :-1], dtype=torch.float32)
    else:
        return (torch.tensor(df.values[:, :-1], dtype=torch.float32),
                torch.tensor(df.values[:, -1], dtype=torch.float32), df)


def get_train_test_data(month_path, train_window=3):
    begin_ind = month_list.index(month_path)
    train_paths = month_list[begin_ind: begin_ind + train_window]
    test_path = month_list[begin_ind + train_window]
    return [read_month_data(p, is_train=True) for p in train_paths], \
           read_month_data(test_path)


def train_model(pair_data):
    model.train()
    data_loader = DataLoader(pair_data, shuffle=True, batch_size=512)
    losses, flag = 0, 0
    for x_i, x_j in tqdm(data_loader):
        sig = model(x_i, x_j)
        l = torch.ones_like(sig)
        loss = loss_func(sig, l)
        opt.zero_grad()
        loss.backward()
        opt.step()
        losses += loss.item()
        flag += 1
    return losses / flag


def val_model(pair_data):
    model.eval()
    data_loader = DataLoader(pair_data, shuffle=True, batch_size=512)
    losses = 0
    for x_i, x_j in tqdm(data_loader):
        sig = model(x_i, x_j)
        # print(sig.view(-1))
        losses += (sig.view(-1) > 0.5).sum().item()
    return losses / len(pair_data)


def eval_model(val_data):  # top stk - bot stk
    model.eval()
    pred_ = model.model(val_data[0]).view(-1)
    sort_p = val_data[1][pred_.argsort()]
    return (sort_p[-100:] + sort_p[:100]).mean(), pred_.data.numpy()


# %%
train_window = 6
if not os.path.exists('./ranknet_ed2'): os.mkdir('./ranknet_ed2')
for ind in range(len(month_list) - train_window - 1):
    train_data, test_data = get_train_test_data(month_list[ind],
                                                train_window=train_window)
    date = test_path = month_list[ind + train_window][-10:]
    model = RankNet(16, drop_out=0.5)
    loss_func = nn.BCELoss()
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    pair_train_data = pairwise_data(train_data, step=100)
    pair_test_data = pairwise_data([test_data[0]], step=100)

    for i in range(5):
        train_loss = train_model(pair_train_data)
        eval_acc = val_model(pair_test_data)
        # eval_acc = eval_model(val_data)
        ret_test, pred_ = eval_model(test_data)
        print('train loss: %.4f | val acc: %.4f | test return: %.4f' % (train_loss,
                                                                        eval_acc, ret_test))
        test_data[2]['pred_%s' % i] = pred_
    test_data[2].to_csv('./ranknet_ed2/%s' % date)
    print('================')
