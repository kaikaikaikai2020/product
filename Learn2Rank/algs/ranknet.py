import torch
import torch.nn as nn
import numpy as np
from torch.utils import data


# %%
class RankNet(nn.Module):

    def __init__(self, num_feature, hidden_size=64, drop_out=0.5):
        super(RankNet, self).__init__()
        hidden_2 = int(np.sqrt(hidden_size))

        self.model = nn.Sequential(
            nn.Linear(num_feature, hidden_size),
            nn.Dropout(drop_out),
            nn.ReLU(),
            nn.Linear(hidden_size, hidden_2),
            nn.Dropout(drop_out),
            nn.ReLU(),
            nn.Linear(hidden_2, 1),
        )
        self.output_sig = nn.Sigmoid()

    def forward(self, input_1, input_2):
        s1 = self.model(input_1)
        s2 = self.model(input_2)
        out = self.output_sig(s1 - s2)

        return out


# %%
class pairwise_data(data.Dataset):

    def __init__(self, order_xs, step=50):
        self.xs = order_xs
        self.x = torch.cat(order_xs, axis=0)
        self.step = step

        self.lens = []
        flag = 0
        for l in self.xs:
            self.lens.append(flag)
            flag += l.shape[0]

        self.inds = []
        for (order_x, base_len) in zip(self.xs, self.lens):
            inds = [i for i in self.pairwise_ind_iter(order_x, base_len)]
            self.inds.extend(inds)

    def pairwise_ind_iter(self, order_x, base_len):  # generate pair data
        len_n = order_x.shape[0]
        for i in range(len_n - self.step):
            for j in range(i + self.step, len_n):
                if i < 100 or j > len_n - 100:
                    yield (i + base_len, j + base_len)

    def __getitem__(self, item):
        high, low = self.inds[item]
        return self.x[high], self.x[low]

    def __len__(self):
        return len(self.inds)


# %%
if __name__ == '__main__':

    n_feature = 16
    model = RankNet(16)
    x_1 = torch.randn(size=(1000, 16))
    x_2 = torch.randn(size=(1000, 16))

    score = model(x_1, x_2)
    train_data = pairwise_data([x_1, x_2, x_2])
    data_loader = data.DataLoader(train_data, shuffle=True, batch_size=256)
    for i in data_loader:
        print(i)
