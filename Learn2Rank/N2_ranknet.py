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
# %%
month_list = glob('L:\Dropbox\Dropbox\project folder from my asua computer\Project\Adaboost\month_data\*.csv')
month_list.sort()


# %%
def read_month_data(path, is_train=False):
    df = pd.read_csv(path, index_col=0)
    df.sort_values('real_return', ascending=False, inplace=True)
    df_ = df.iloc[:, 3:-3]
    if is_train:
        return torch.tensor(df_.values[:, :-1], dtype=torch.float32)
    else:
        return (torch.tensor(df_.values[:, :-1], dtype=torch.float32),
                torch.tensor(df_.values[:, -1], dtype=torch.float32),
                df.loc[:, ['ticker', 'tradeDate', 'real_return', 'market_value']])


def get_train_test_data(month_path, train_window=3):
    begin_ind = month_list.index(month_path)
    train_paths = month_list[begin_ind: begin_ind + train_window]
    val_path = month_list[begin_ind + train_window]
    test_path = month_list[begin_ind + train_window + 1]
    return [read_month_data(p, is_train=True) for p in train_paths], \
           read_month_data(val_path), read_month_data(test_path)


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
    return (sort_p[-100:] - sort_p[:100]).mean(), pred_.data.numpy()


# %%
train_window = 6
n_round = 10
patience = 3
if not os.path.exists('./ranknet'): os.mkdir('./ranknet')
for ind in range(len(month_list) - train_window-1):
    print(ind)
    train_data, val_data, test_data = get_train_test_data(month_list[ind],
                                                          train_window=train_window)
    date = month_list[ind + train_window + 1][-11:]
    model = RankNet(244, drop_out=0.5)
    loss_func = nn.BCELoss()
    opt = torch.optim.Adam(model.parameters(), lr=1e-5)
    pair_train_data = pairwise_data(train_data, step=100)
    pair_val_data = pairwise_data([val_data[0]], step=100)

    eval_accs = []
    acc = 0
    for i in range(n_round):
        train_loss = train_model(pair_train_data)
        eval_acc = val_model(pair_val_data)
        # test_ret, _ = eval_model(test_data)
        # eval_acc = eval_model(val_data)
        test_ret, pred_ = eval_model(test_data)
        print('train loss: %.4f | val acc: %.4f | test return: %.4f' % (train_loss,
                                                                        eval_acc, test_ret))
        eval_accs.append((eval_acc, pred_))
        if eval_acc > acc:
            acc = eval_acc
            flag = 0
        flag += 1
        if flag > patience:
            break

    test_data[2]['pred'] = eval_accs[np.array([i[0] for i in eval_accs]).argmax()][1]  # get the best prediction
    test_data[2].to_csv('./ranknet/%s' % date)
    print('================%s' % date)
