import pandas as pd
import numpy as np
import torch
from torch import nn
from torch.utils.data import TensorDataset, DataLoader
from algs.listnet import *
from glob import glob
from tqdm import tqdm
import copy

torch.random.manual_seed(2021)

# %%
month_list = glob('/home/dinker/桌面/codes/Adaboost/month_data/*.csv')
month_list.sort()


# %%
def read_month_data(path):
    df = pd.read_csv(path, index_col=0)
    return (torch.tensor(df.values[:, :-1], dtype=torch.float32),
            torch.tensor(df.values[:, -1], dtype=torch.float32))


def get_train_val_test_data(month_path, train_window=3):
    begin_ind = month_list.index(month_path)
    train_paths = month_list[begin_ind: begin_ind + train_window]
    val_path = month_list[begin_ind + train_window]
    test_path = month_list[begin_ind + train_window + 1]
    return [read_month_data(p) for p in train_paths], \
           read_month_data(val_path), \
           read_month_data(test_path)


# %%
r_ind = np.random.randint(len(month_list))
print(r_ind)
r_ind = 17
train_data, val_data, test_data = get_train_val_test_data(month_list[r_ind], train_window=6)

# %%
model = ListNet(16, drop_out=0.5)
loss_func = ListMLE_loss()
opt = torch.optim.Adam(model.parameters(), lr=1e-3)

train_td = TensorDataset(torch.cat([i[0] for i in train_data], axis=0),
                         torch.cat([i[1] for i in train_data], axis=0))


# pair_train_data = pairwise_data(train_data, step=200)
# pair_test_data = pairwise_data([val_data[0]], step=200)
# train_size = int(0.8 * len(pair_data))
# test_size = len(pair_data) - train_size
# train_pair, test_pair = random_split(pair_data, (train_size, test_size))

# %%

def train_model(train_td):
    model.train()
    data_loader = DataLoader(train_td, shuffle=True, batch_size=256)
    losses, flag = 0, 0
    for x_i, y_i in tqdm(data_loader):
        sig = model(x_i)
        loss = loss_func(sig, y_i)
        opt.zero_grad()
        loss.backward()
        opt.step()
        losses += loss.item()
        flag += 1
    return losses / flag


def val_model(val_data):
    model.eval()
    sig = model(val_data[0])
    score = ndcg(sig, val_data[1])
    return score


def eval_model(test_data, model):  # top stk - bot stk
    model.eval()
    pred_ = model(test_data[0]).view(-1)
    sort_p = test_data[1][pred_.argsort()]
    return (sort_p[-100:] - sort_p[:100]).mean()


# %%
flag = True
for i in range(100):
    train_loss = train_model(train_td)
    eval_acc = val_model(val_data)
    ret_test = eval_model(test_data, model)
    print('train loss: %.4f | val ndcg@100: %.4f | test return: %.4f' % (train_loss,
                                                                         eval_acc, ret_test))
    # if flag:
    #     best_model = copy.deepcopy(model)
    #     flag = False
    # best_score = eval_model(test_data, best_model)
    # print(best_score)


# %%
# model.eval()
# data_loader = DataLoader(test_pair, shuffle=True, batch_size=512)
# losses = 0
# for x_i, x_j in tqdm(data_loader):
#     sig = model(x_i, x_j)
#     print(sig.view(-1))
#     losses += (sig.view(-1) > 0.5).sum()
