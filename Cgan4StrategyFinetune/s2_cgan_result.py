import torch
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# %%
g_model = torch.load('./g_model.pkl')
data = pd.read_csv('./data/hs_300.csv')

# %%
s = data['CHGPct'].values
window_size, n = 252, data.shape[0]
c = int((n - window_size) * .8)
x_data = np.array([s[i: i + window_size] for i in range(n - window_size)])
y_data = np.array([s[i + window_size] for i in range(n - window_size)]).reshape((-1, 1))


# %% simulate
def cgan_sim(t='train'):
    if t == 'train':
        x = torch.tensor(x_data[:c], dtype=torch.float32)
        z = torch.randn(size=(c, 252))
        return torch.cat([z, x], axis=1)
    else:
        x = torch.tensor(x_data[c:], dtype=torch.float32)
        b, p = x.shape
        z = torch.randn(size=(b, p))
        return torch.cat([z, x], axis=1)


# %%
y_gan = torch.cat([g_model(cgan_sim('train')).data for i in range(300)], axis=1)

# %%
plt.figure(figsize=(16, 8))
plt.plot(y_data[:c], linewidth=0.2)
plt.plot(y_gan, linewidth=0.1)
plt.show()

# %%
plt.figure(figsize=(16, 8))
plt.plot((y_data[:c] + 1).cumprod(axis=0), linewidth=0.3)
plt.plot((y_gan + 1).cumprod(axis=0), linewidth=0.1)
plt.show()

# %%
# y_gan = torch.cat([g_model(cgan_sim('test')).data for i in range(300)], axis = 1)
#
# plt.plot((y_data[c:]+ 1).cumprod(axis = 0), linewidth=0.3)
# plt.plot((y_gan + 1).cumprod(axis = 0), linewidth=0.1)
# plt.show()
#
# # %%
# p = y_data[c:]
# plt.plot((p[y_gan.mean(axis = 1) > 0] + 1).cumprod())
# plt.show()
# %%
# plt.plot((y_data[:c] + 1).cumprod(axis = 0), linewidth=0.3)
# plt.plot((y_gan + 1).cumprod(axis = 0), linewidth=0.1)
# plt.show()
