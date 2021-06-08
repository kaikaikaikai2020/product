import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt


# %%
class bubble_boost:

    def __init__(self, y, boosting_round=500):
        delta_y = y[1:] - y[:-1]
        n = len(y)
        self.delta_y, self.delta_ym1 = delta_y[1:], delta_y[:-1]
        self.y, self.ym1, self.y1 = y[2:], y[1:-1], y[1]
        self.tau = int((0.01 + 1.8 / np.sqrt(n)) * n)
        self.boosting_round = boosting_round

    def get_delta_info(self, y): # use to get base adf
        delta_y = y[1:] - y[:-1]
        delta_y, delta_ym1 = delta_y[1:], delta_y[:-1]
        y, ym1, y1 = y[2:], y[1:-1], y[1]
        return delta_y, ym1, delta_ym1

    def get_adf(self, delta_y, ym1, delta_ym1):
        x = sm.add_constant(np.vstack([ym1, delta_ym1]).T)
        result = sm.OLS(delta_y, x).fit()
        return result.tvalues

    def get_max_adf(self, delta_y, ym1, delta_ym1):
        n = len(delta_y)
        adfs = np.zeros(shape=(n - self.tau,))
        for i in range(n - self.tau):
            result = self.get_adf(delta_y[i:], ym1[i:], delta_ym1[i:])
            adfs[i] = result[1]
        return adfs.max()

    def boosting(self):
        x = sm.add_constant(self.delta_ym1)
        result = sm.OLS(self.delta_y, x).fit()
        eta = result.resid

        # boosting
        boosting_result = []
        for _ in range(self.boosting_round):
            w = np.random.normal(0, 1, size=(len(eta),))
            boost_ind = np.sort(np.random.randint(0, len(eta), size=(len(eta))))
            boost_delta_y = sm.add_constant(self.delta_ym1[boost_ind]) @ result.params + eta[boost_ind] * w
            # new_y = self.y1 + np.cumsum(boost_delta_y)
            new_delta_y = boost_delta_y
            new_ym1, new_delta_ym1 = self.ym1[boost_ind], self.delta_ym1[boost_ind]
            boosting_result.append([new_delta_y, new_ym1, new_delta_ym1])
        return boosting_result

    def fit(self):
        true_madf = self.get_max_adf(self.delta_y, self.ym1, self.delta_ym1)
        print('true madf has done')
        boosting_result = self.boosting()
        print('boosting has done')
        max_adfs = np.zeros(shape=(self.boosting_round))
        for ind, (delta_y, ym1, delta_ym1) in enumerate(boosting_result):
            max_adf = self.get_max_adf(delta_y, ym1, delta_ym1)
            max_adfs[ind] = max_adf
            if ind % 20 == 0:
                print('boosting %s max adf has done' % ind)
        return true_madf, max_adfs

    def generate_tnesor_xy(self):
        boosting_result = self.boosting()
        y_tensor = np.zeros(shape=(self.boosting_round, len(self.y), 1))
        x_tensor = np.zeros(shape=(self.boosting_round, len(self.y), 3))
        for ind, (delta_y, ym1, delta_ym1) in enumerate(boosting_result):
            y_tensor[ind] = delta_y.reshape(-1, 1)
            x = sm.add_constant(np.vstack([ym1, delta_ym1]).T)
            x_tensor[ind] = x
        return x_tensor, y_tensor

    def get_tensor_tratio(self, x_tensor, y_tensor):
        n = y_tensor.shape[1]
        x_trans = x_tensor.transpose(0, 2, 1)
        xTx = np.matmul(x_trans, x_tensor)
        xTx_inv = np.linalg.inv(xTx)
        params = np.matmul(np.matmul(xTx_inv, x_trans), y_tensor)
        e_tensor = (y_tensor - np.matmul(x_tensor, params)).reshape(self.boosting_round, -1)
        t_vals = params[:, 1, :].ravel() / np.sqrt(
            (xTx_inv[:, 1, 1].ravel() * (e_tensor * e_tensor).sum(axis=1) / (n - 3)))
        return t_vals

    def fast_fit(self):
        # true_madf = self.get_max_adf(self.delta_y, self.ym1, self.delta_ym1)
        # print('true madf has done')
        x_tensor, y_tensor = self.generate_tnesor_xy()
        # print('boosting tensor has done')
        n = len(self.delta_y)
        adfs = np.zeros(shape=(n - self.tau, self.boosting_round))
        for i in range(n - self.tau):
            x_tensor_ = x_tensor[:, i:, :]
            y_tensor_ = y_tensor[:, i:, :]
            t_vals = self.get_tensor_tratio(x_tensor_, y_tensor_)
            adfs[i] = t_vals
        # print('boosting adfs has done')
        return adfs.max(axis=0)

    def get_true_adf(self):
        true_madf = self.get_max_adf(self.delta_y, self.ym1, self.delta_ym1)
        return true_madf



# %%
if __name__ == '__main__':
    df = pd.read_csv('./data/indu_week.csv', index_col=0)
    indu_id = df['industryID1'].unique().tolist()

    # %%
    sub_df = df[df['industryID1'] == indu_id[0]].copy()
    sub_df.index = sub_df.pop('endDate')
    sub_df.sort_index(inplace=True)
    sub_df = sub_df.iloc[4:]  # 对刚上的指数进行一些基本操作

    # %%
    y, delta_y = sub_df['closePrice'].values, sub_df['chg'].values
    model = bubble_boost(y)
    model.generate_tnesor_xy()
    adfs = model.fast_fit()
    # a, b = model.fit()
    #
    # %%
    # x = sm.add_constant(model.delta_ym1)
    # result = sm.OLS(model.delta_y, x).fit()
    # eta = result.resid
    # w = np.random.normal(0, 1, size=(len(eta),))
    # pred_delta = x @ result.params + eta * w
    #
    # # %%
    # plt.plot(eta)
    # plt.show()
    #
    #
    # # %%
    # plt.plot(model.y1 + model.delta_y.cumsum())
    # plt.plot(model.y1 + pred_delta.cumsum())
    # plt.show()
