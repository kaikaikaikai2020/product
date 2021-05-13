import torch
import torch.nn as nn
from tqdm import tqdm

# %%
class Generator(nn.Module):

    def __init__(self, input_dim, output_dim, hidden_dim=100):
        super().__init__()
        self.main_model = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, output_dim)
        )

    def forward(self, x):
        g = self.main_model(x)
        return g


class Discriminator(nn.Module):

    def __init__(self, input_dim, output_dim, hidden_dim=100):
        super().__init__()
        self.main_model = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, output_dim),
            nn.Sigmoid()
        )

    def forward(self, x):
        d = self.main_model(x)
        return d


class CGAN_model:

    def __init__(self, g_model, d_model, epochs=20000, opt_lr=1e-3, fake_size=252):
        self.g_model = g_model
        self.d_model = d_model
        self.opt_g = torch.optim.Adam(self.g_model.parameters(), lr=opt_lr)
        self.opt_d = torch.optim.Adam(self.d_model.parameters(), lr=opt_lr)
        self.loss_func = nn.BCELoss()
        self.fake_size = fake_size
        self.epochs = epochs

    # Update D network: maximize log(D(x)) + log(1 - D(G(z)))
    def update_Dnet(self, real_data, fake_data):
        self.opt_d.zero_grad()
        y_real = self.d_model(real_data)
        err_real = self.loss_func(y_real, torch.ones_like(y_real))
        y_fake = self.d_model(fake_data)
        err_fake = self.loss_func(y_fake, torch.zeros_like(y_real))
        err_d = err_fake + err_real
        # print(err_real.item(), err_fake.item())
        err_d.backward()
        self.opt_d.step()

    # Update G network: maximize log(D(G(z)))
    def update_Gnet(self, fake_g, condition_x):
        self.opt_g.zero_grad()
        y_gen = self.d_model(torch.cat([fake_g, condition_x], axis=1))
        err_g = self.loss_func(y_gen, torch.ones_like(y_gen))
        err_g.backward()
        self.opt_g.step()

    def generate_data(self, condition_x):
        b, n = condition_x.shape
        z = torch.randn(size=(b, self.fake_size))
        return torch.cat([z, condition_x], axis=1)

    def eval_result(self, real_y, condition_x):
        fake_y = self.g_model(self.generate_data(condition_x))
        score = torch.sum((real_y - fake_y) ** 2)
        return score.item()

    def train_step(self, condition_x, real_y):
        real_data = torch.cat([real_y, condition_x], axis=1)
        gen_data = self.generate_data(condition_x)

        fake_g = self.g_model(gen_data)
        fake_data = torch.cat([fake_g, condition_x], axis=1)
        self.update_Dnet(real_data, fake_data)

        # last time the gradient info has used
        fake_g = self.g_model(gen_data)
        self.update_Gnet(fake_g, condition_x)

    def main_loop(self, dataloader):
        scores = []
        for _ in tqdm(range(self.epochs)):
            score = 0
            n = 0
            for i, data in enumerate(dataloader):
                condition_x, real_y = data
                self.train_step(condition_x, real_y)
                score += self.eval_result(real_y, condition_x)
                n += real_y.shape[0]
            scores.append(score / len(dataloader))
            print(score / n)
        return scores


# %%
if __name__ == '__main__':
    x_size, p_size = 252, 252
    g_model = Generator(x_size + p_size, 1)
    d_model = Discriminator(x_size + 1, 1)
    cgan = CGAN_model(g_model, d_model, fake_size=252)
    x = torch.randn(size=(256, 252))
    y = torch.randn(size=(256, 1))
    cgan.train_step(x, y)
