import torch
import torch.nn as nn
import numpy as np
from torch.utils import data


# %%
class ListNet(nn.Module):

    def __init__(self, num_feature, hidden_size=64, drop_out=0.5):
        super(ListNet, self).__init__()
        hidden_2 = int(np.sqrt(hidden_size))

        self.model = nn.Sequential(
            nn.Linear(num_feature, hidden_size),
            nn.Dropout(drop_out),
            nn.ReLU(),
            nn.Linear(hidden_size, hidden_2),
            nn.Dropout(drop_out),
            nn.ReLU(),
            nn.Linear(hidden_2, 1)
        )

    def forward(self, input_x):
        out = self.model(input_x)
        return out


class ListMLE_loss(nn.Module):

    def __init__(self, eps=1e-5):
        super(ListMLE_loss, self).__init__()
        self.eps = eps

    def forward(self, y_pred, y_true):
        y_sort, indices = y_true.sort(descending=True, dim=0)
        indices = indices.view(-1)
        y_sort_pred = y_pred[indices]
        cumloss = torch.cumsum(y_sort_pred.exp().flip(dims=[0]), dim=0).flip(dims=[0])
        observed_loss = torch.log(cumloss + self.eps) - y_sort_pred
        return observed_loss.sum()


def ndcg(ys_true, ys_pred, n=100):
    # quantile
    ys_true, ys_pred = ys_true.data.numpy().ravel(), ys_pred.data.numpy().ravel()
    ys_true[ys_true.argsort()] = np.arange(1, ys_true.shape[0] + 1) / ys_true.shape[0]
    ys_pred[ys_pred.argsort()] = np.arange(1, ys_pred.shape[0] + 1) / ys_pred.shape[0]
    # assert (ys_true > 0).all()

    def dcg(ys_true, ys_pred):
        argsort = ys_pred.argsort()[::-1]
        ys_true_sorted = ys_true[argsort]
        ret = 0
        for i, l in enumerate(ys_true_sorted[:n], 1):
            ret += (2 ** l - 1) / np.log2(1 + i)
        return ret

    ideal_dcg = dcg(ys_true, ys_true)
    pred_dcg = dcg(ys_true, ys_pred)
    return (pred_dcg / ideal_dcg)


# %%
if __name__ == '__main__':
    ws = torch.randn(16, 1)
    x_train = torch.randn(1000, 16)
    x_val = torch.randn(1000, 16)
    y_train = torch.mm(x_train, ws)
    y_val = torch.mm(x_val, ws)

    #
    model = ListNet(16)
    opt = torch.optim.Adam(model.parameters(), lr=0.01)
    loss_func = ListMLE_loss()
    for _ in range(10):
        model.train()
        y_pred = model(x_train)
        loss = loss_func(y_pred, y_train)
        opt.zero_grad()
        loss.backward()
        opt.step()
        model.eval()
        print('train loss: %.4f| val ndcg@100: %.4f' % (loss.item(), ndcg(model(x_val), y_val)))

