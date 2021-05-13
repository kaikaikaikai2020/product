import torch
import torch.nn as nn

# %%
model_1 = nn.Linear(10, 10)
model_2 = nn.Linear(10, 1)

opt = torch.optim.Adam(model_2.parameters())
loss_func = nn.MSELoss()
x1_w = model_1.weight
# %%
x = torch.randn(size=(265, 10))
y = torch.randn(size=(265, 1))
for i in range(2):
    x_1 = model_1(x)
    x_2 = model_2(x_1)
    loss = loss_func(x_2, y)
    assert torch.all(model_1.weight == x1_w)
    opt.zero_grad()
    loss.backward()
    opt.step()