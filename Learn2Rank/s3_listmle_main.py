import pandas as pd
import os
import torch
from torch import nn
from torch.utils.data import TensorDataset, DataLoader
from algs.listnet import *
from glob import glob
from tqdm import tqdm
import copy

torch.random.manual_seed(2021)

# %%

month_list = glob('./month/*.csv')
month_list.sort()


# %%
def read_month_data(path):
    df = pd.read_csv(path, dtype={'ts_code': str}, index_col='ts_code')
    return (torch.tensor(df.values[:, :-1], dtype=torch.float32),
            torch.tensor(df.values[:, -1], dtype=torch.float32), df)


def get_train_val_test_data(month_path, train_window=3):
    begin_ind = month_list.index(month_path)
    train_paths = month_list[begin_ind: begin_ind + train_window]
    val_path = month_list[begin_ind + train_window]
    test_path = month_list[begin_ind + train_window + 1]
    return [read_month_data(p) for p in train_paths], \
           read_month_data(val_path), \
           read_month_data(test_path)


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
    return (sort_p[-100:] - sort_p[:100]).mean(), pred_.data.numpy()


# %%
train_window = 6
patience = 10
if not os.path.exists('./listmle_ed2'): os.mkdir('./listmle_ed2')
for ind in range(len(month_list) - train_window - 1):
    date = test_path = month_list[ind + train_window][-10:]
    train_data, val_data, test_data = get_train_val_test_data(month_list[ind], train_window=train_window)
    model = ListNet(16, drop_out=0.5)
    loss_func = ListMLE_loss()
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    train_td = TensorDataset(torch.cat([i[0] for i in train_data], axis=0),
                             torch.cat([i[1] for i in train_data], axis=0))

    max_ndcg = -np.inf
    p_flag = 0
    for i in range(100):
        train_loss = train_model(train_td)
        eval_acc = val_model(val_data)
        ret_test, _ = eval_model(test_data, model)
        print('train loss: %.4f | val ndcg@100: %.4f | test return: %.4f' % (train_loss,
                                                                             eval_acc, ret_test))
        if eval_acc > max_ndcg:
            max_ndcg = eval_acc
            best_model = copy.deepcopy(model)
            p_flag = 0
        p_flag += 1
        if p_flag > patience:
            break
    _, pred_ = eval_model(test_data, best_model)
    test_data[2]['pred'] = pred_
    test_data[2].to_csv('./listmle_ed2/%s' % date)
    print('========== %s has done' % date)

# %%
# model.eval()
# data_loader = DataLoader(test_pair, shuffle=True, batch_size=512)
# losses = 0
# for x_i, x_j in tqdm(data_loader):
#     sig = model(x_i, x_j)
#     print(sig.view(-1))
#     losses += (sig.view(-1) > 0.5).sum()
