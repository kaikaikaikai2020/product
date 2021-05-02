import numpy as np



# %%
class macd_direction:

    def __init__(self, fast=12, slow=26, trigger=9, theta=0.005):
        self.fast = fast
        self.slow = slow
        self.trgger = trigger
        self.theta = theta

    def macd_(self, s):
        s_fast = s.ewm(adjust=False, alpha=2 / (self.fast + 1), ignore_na=True).mean()
        s_slow = s.ewm(adjust=False, alpha=2 / (self.slow + 1), ignore_na=True).mean()
        s_dif = s_fast - s_slow
        s_trig = s_dif.ewm(adjust=False, alpha=2 / (self.trgger + 1), ignore_na=True).mean()
        return s_dif, s_trig

    def integral(self, s):
        s_dif, signal = self.macd_(s)
        s_macd = (s_dif - signal).values
        inter_arr = np.zeros_like(s_dif)
        val = 0
        for ind, i in enumerate(s_macd):
            if np.sign(val) == 0:
                val += i
            if np.sign(val) == np.sign(i):
                val += i
            else:
                val = 0
            inter_arr[ind] = val
        return inter_arr

    def integal_macd(self, s):
        inter = self.integral(s)
        thread = s * self.theta
        macd_dir = np.zeros_like(s)
        temp_flag = 1
        for ind, (i, t) in enumerate(zip(inter, thread)):
            if i > t or (temp_flag == 1 and i > -t):
                inter_dir = 1
            elif i < -t or (temp_flag == -1 and i < t):
                inter_dir = -1
            temp_flag = inter_dir
            macd_dir[ind] = inter_dir
        return macd_dir


# %%
if __name__ == '__main__':
    import pandas as pd
    import matplotlib.pyplot as plt

    df = pd.read_csv('./data/hs_300.csv', index_col=0)
    s = df['closeIndex']
    macd_dir = macd_direction().integal_macd(s)

    # plot
    plt.figure(12, figsize=(16, 10))
    plt.subplot(211)
    plt.plot(s[1500:2000])
    plt.subplot(212)
    plt.plot(macd_dir[1500:2000])
    plt.show()

# %%
# s_slow = s.ewm(adjust=False, alpha=2 / (26 + 1), ignore_na=True).mean()
# s_fast = s.ewm(adjust=False, alpha=2 / (12 + 1), ignore_na=True).mean()
# s_dif = s_fast - s_slow
# s_trig = s_dif.ewm(adjust=False, alpha=2 / (9 + 1), ignore_na=True).mean()
# s_macd = s_dif - s_trig
#
# # %%
# val = 0
# for i in s_macd:
#     if np.sign(val) == 0:
#         val += i
#     if np.sign(val) == np.sign(i):
#         val += i
#     else:
#         val = 0
#
#
# # %%
# plt.figure(12, figsize=(16, 10))
# plt.subplot(211)
# plt.plot(s[350:500])
# plt.plot(s_fast[350:500])
# plt.plot(s_slow[350:500])
# plt.subplot(212)
# plt.plot(s_dif[350:500])
# plt.plot(s_trig[350:500])
# # plt.plot(s_macd[350:500])
# plt.show()
#
# # %%
# s.ewm(adjust=False, alpha=2 / (12 + 1), ignore_na=True).mean()
