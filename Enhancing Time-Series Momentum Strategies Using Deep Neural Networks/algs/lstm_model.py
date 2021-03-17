import torch
from torch import nn
import numpy as np


# %% lstm model
class lstm_model(nn.Module):

    def __init__(self, input_feat=20, hidden_dim=5, dropout=0.3):
        super().__init__()
        self.dropout = dropout
        self.lstm = nn.LSTM(input_size=input_feat,
                            hidden_size=hidden_dim,
                            batch_first=True,
                            num_layers=1)
        self.dropout_layer = nn.Dropout(self.dropout)
        self.linear = nn.Linear(in_features=hidden_dim, out_features=1)

    def forward(self, x):
        lstm_st, (c_n, h_n) = self.lstm(x)
        lstm_st, _, _ = self.dropout_layer(lstm_st), \
                        self.dropout_layer(c_n), self.dropout_layer(h_n),
        output = self.linear(lstm_st[:, -1, :])
        return output


class avg_return_loss(nn.Module):

    def __init__(self):
        super().__init__()
        # trans position to (-1, 1)
        self.position_trans = nn.Tanh()

    def forward(self, x, y):
        pos = self.position_trans(x)
        return - (pos * y).mean()


class reg_loss(nn.Module):

    def __init__(self):
        super().__init__()
        self.mse = nn.MSELoss(reduction='mean')

    def forward(self, x, y):
        return self.mse(x, y)


class binary_loss(nn.Module):

    def __init__(self):
        super().__init__()
        # trans x to (0, 1)
        self.prob = nn.Sigmoid()
        # trans y to {0, 1}
        self.label_fn = lambda x: torch.where(x >= 0, torch.ones_like(x, dtype=torch.float),
                                              torch.zeros_like(x, dtype=torch.float))
        # define loss func
        self.loss = nn.BCELoss()

    def forward(self, x, y):
        pr_x = self.prob(x)
        pr_y = self.label_fn(y)
        return self.loss(pr_x, pr_y)


class sharpe_loss(nn.Module):

    def __init__(self):
        super().__init__()
        self.position = nn.Tanh()

    def forward(self, x, y):
        pos = self.position(x)
        rt_ser = pos * y
        shp = rt_ser.sum() / ((rt_ser ** 2).sum() - rt_ser.sum() ** 2)
        return - shp

class EarlyStopping:
    """Early stops the training if validation loss doesn't improve after a given patience."""
    def __init__(self, patience=7, verbose=False, delta=0, checkpoint = 'checkpoint.pt'):
        """
        Args:
            patience (int): How long to wait after last time validation loss improved.
                            Default: 7
            verbose (bool): If True, prints a message for each validation loss improvement.
                            Default: False
            delta (float): Minimum change in the monitored quantity to qualify as an improvement.
                            Default: 0
        """
        self.patience = patience
        self.verbose = verbose
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.val_loss_min = np.Inf
        self.delta = delta
        self.checkpoint = checkpoint

    def __call__(self, val_loss, model):

        score = -val_loss

        if self.best_score is None:
            self.best_score = score
            self.save_checkpoint(val_loss, model)
        elif score < self.best_score + self.delta:
            self.counter += 1
            # print(f'EarlyStopping counter: {self.counter} out of {self.patience}')
            if self.counter >= self.patience:
                self.early_stop = True
        else:
            self.best_score = score
            self.save_checkpoint(val_loss, model)
            self.counter = 0

    def save_checkpoint(self, val_loss, model):
        '''Saves model when validation loss decrease.'''
        if self.verbose:
            print(f'Validation loss decreased ({self.val_loss_min:.6f} --> {val_loss:.6f}).  Saving model ...')
        torch.save(model.state_dict(), self.checkpoint)
        self.val_loss_min = val_loss

# %%
if __name__ == '__main__':
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    x = torch.rand(size=(256, 50, 20))
    x = x.to(device)
    y = torch.rand(size=(256, 1))
    y = y.to(device)

    model = lstm_model()
    model.to(device)
    pred_y = model(x)
    for loss_func in [reg_loss(), binary_loss(), avg_return_loss(), sharpe_loss()]:
        loss_func = loss_func.to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=5e-4)
        pred_y = model(x)
        loss = loss_func(pred_y, y)
        loss.backward()
        optimizer.step()

        print(loss)
        print('test is pass')
