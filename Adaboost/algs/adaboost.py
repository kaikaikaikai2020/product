import numpy as np
import pandas as pd


# %%
class adaboost_weak_cls():

    def __init__(self, n_models=30, Q=5, w=None):
        self.Q = Q
        self.w = w
        self.n_models = n_models
        self.pred_list = []  # 用来做预测

    def get_group(self, X):
        x_group = np.apply_along_axis(lambda x:
                                      pd.qcut(x, self.Q, range(self.Q)), 0, X)
        return x_group

    def get_error(self, x, y, w):  # group_x
        y = y.ravel()
        sig_pos = y == 1
        sig_neg = y == -1
        wqs_pos = np.zeros(shape=self.Q)
        wqs_neg = np.zeros(shape=self.Q)

        for i in range(self.Q):
            sig_sample = x == i
            wq_pos = w[sig_pos & sig_sample].sum()
            wq_neg = w[sig_neg & sig_sample].sum()
            wqs_pos[i] = wq_pos
            wqs_neg[i] = wq_neg
        return np.sqrt(wqs_pos * wqs_neg).sum()

    def get_confidence(self, x, y, w):
        y = y.ravel()
        sig_pos = y == 1
        sig_neg = y == -1
        hq_update = np.zeros_like(w)
        hq_dict = {}

        for i in range(self.Q):
            sig_sample = x == i
            wq_pos = w[sig_pos & sig_sample].sum()
            wq_neg = w[sig_neg & sig_sample].sum()
            hq = 1 / 2 * np.log(wq_pos / wq_neg)
            hq_update[sig_sample] = hq
            hq_dict[i] = hq

        return hq_update, hq_dict


    def fit(self, X, y):
        n, m = X.shape
        y = y.ravel()
        if not self.w: self.w = np.ones(n) / n
        x_group = self.get_group(X)
        hq_sum = np.zeros_like(self.w)
        for i in range(self.n_models):
            error = np.apply_along_axis(lambda x: self.get_error(x, y=y, w=self.w),
                                        0, x_group)
            index_k = error.argmin()
            hq, hq_dict = self.get_confidence(x_group[:, index_k], y, self.w)
            self.pred_list.append((index_k, hq_dict))
            w_update = self.w * np.exp(- hq * y)
            self.w = w_update / w_update.sum()
            hq_sum += hq

        return hq_sum.reshape(-1, 1)

    def predict(self, X):
        n, m = X.shape
        pred_sum = np.zeros(n)
        x_group = self.get_group(X)
        for s in self.pred_list:
            k, hq_dict = s
            x = x_group[:, k]
            pred_hq = np.zeros(n)
            for i in range(self.Q):
                pred_hq[x == i] = hq_dict[i]
            pred_sum += pred_hq
        return pred_sum.reshape(-1, 1)


# %%
if __name__ == '__main__':
    X = np.random.random(size=(1000, 20))
    y = np.random.randint(0, 2, size=(1000, 1)) * 2 - 1
    model = adaboost_weak_cls()
    pred = model.fit(X, y)
    pred_x = model.predict(X)
