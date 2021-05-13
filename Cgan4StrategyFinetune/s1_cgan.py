import pandas as pd
from algs.cgan import *
from torch.utils.data import DataLoader, Dataset, dataset
import torch

# %%
data = pd.read_csv('./data/hs_300.csv', index_col=0)


class tensor_data(Dataset):

    def __init__(self, x, window_size=252):
        self.x = x
        self.window_size = window_size

    def __getitem__(self, item):
        x_i = torch.FloatTensor(self.x[item: item + self.window_size])
        y_i = torch.FloatTensor([self.x[item + self.window_size]])
        return x_i, y_i

    def __len__(self):
        return len(self.x) - self.window_size


# %%
s = data['CHGPct'].values
# s_std = (s - s.mean()) / s.std()
# s_std = s * 20
print(s.var())
t_data = tensor_data(s)
train_data = dataset.Subset(t_data, range(int(0.8 * len(t_data))))
data_load = DataLoader(train_data, batch_size=256, shuffle=True)

# %%
x_size, p_size = 252, 252
g_model = Generator(x_size + p_size, 1, hidden_dim=50)
d_model = Discriminator(x_size + 1, 1, hidden_dim=50)
cgan = CGAN_model(g_model, d_model, fake_size=252, epochs=10000)
x = torch.randn(size=(256, 252))
y = torch.randn(size=(256, 1))
cgan.main_loop(data_load)

# %%
torch.save(cgan.g_model, './g_model.pkl')
torch.save(cgan.d_model, './d_model.pkl')
