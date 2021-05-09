import numpy as np
from alg_macd_direction import macd_direction
import pandas as pd

# %%
df = pd.read_csv('./data/hs_300.csv', index_col=0)
s = df['closeIndex'].copy()
macd_dir = macd_direction().integal_macd(s).copy()
break_flag = True

# %% macd的区间划分不存在未来信息

while break_flag:
    # print('ok', macd_dir.shape[0], s.shape[0])

    max_dir, min_dir = np.zeros_like(macd_dir), np.zeros_like(macd_dir)
    max_last, min_last = np.zeros_like(macd_dir), np.zeros_like(macd_dir)
    # init
    confirm_up, confirm_down, exp_flag = False, False, False
    b_flag, c_flag = False, False
    b_ind = 0
    t_max, t_min, t_dir = 1, 10e6, macd_dir[0]
    l_max, l_min = s[0], s[0]
    for ind, (s_, dir_) in enumerate(zip(s.values, macd_dir)):
        if t_dir != dir_:
            l_max, l_min = t_max, t_min
            t_max, t_min = s_, s_
            l_ind = b_ind
            b_ind = ind

        # update min max
        if t_max < s_:
            t_max = s_
            # update before max info
            max_dir[b_ind:ind + 1] = t_max
        if t_min > s_:
            t_min = s_
            # update before min info
            min_dir[b_ind:ind + 1] = t_min

        # update exp
        if dir_ == 1 and t_min < l_min and t_max < l_max:  # up-dir exp confirm
            confirm_up = True
        if dir_ == -1 and t_min > l_min and t_max > l_max:
            confirm_down = True
        # if have confirm change dir
        if confirm_up or confirm_down:
            b_flag = True
        if (confirm_up and t_max > l_max) or (confirm_down and t_min < l_min):
            c_flag = True

        if t_dir != dir_ and b_flag and not c_flag:
            exp_flag = True
        if exp_flag:
            macd_dir[l_ind:b_ind] = -macd_dir[l_ind:b_ind]
            print(ind)
            break

        # update max min lmax, lmin
        max_dir[ind] = t_max
        min_dir[ind] = t_min
        max_last[ind] = l_max
        min_last[ind] = l_min

        # 保留上一部dir
        t_dir = dir_

        # all exp has refine
        if ind == macd_dir.shape[0] - 1:
            break_flag = False

# %%
s_df = df[['closeIndex']].copy()
s_df['dir'] = macd_dir
s_df['max_dir'] = max_dir
s_df['min_dir'] = min_dir
s_df['max_last'] = max_last
s_df['min_last'] = min_last
s_minmax = s_df.apply(lambda x: x['max_dir'] if x['dir'] == 1 else x['min_dir'], axis=1)
# %%
import matplotlib.pyplot as plt

plt.figure(figsize=(16, 5))
plt.plot(s[1500:2000])
plt.plot(s_minmax[1500:2000])
plt.show()

# %% 策略exp  closeIndex > max_last表示方向确认，此时 exp 不起作用。没有涉及未来信息。
sign_long = s_df.apply(lambda x: 1 if x['dir'] == 1 and x['closeIndex'] > x['max_last'] else 0, axis=1)
s_df['chg'] = s / s.shift(1) - 1
plt.figure(figsize=(16, 8))
(s_df['chg'] + 1).cumprod().plot(label='index')
(s_df['chg'] * sign_long.shift(1) + 1).cumprod().plot(label='stra1')
# plt.show()

# %% 策略1 在第一次趋势确认时，进行trigger操作。
sig = 0
sigs = np.zeros_like(macd_dir)
s_, t_dir = s[0], macd_dir[0]
ll_max, lll_max, llll_max = max_last[0], max_last[0], max_last[0]
ll_min, lll_min, llll_min = min_last[0], min_last[0], min_last[0]
for ind, (s_, dir_, l_max, l_min) in enumerate(zip(s, macd_dir, max_last, min_last)):

    if dir_ == 1:
        if s_ > ll_max:
            sig = 1
        else:
            sig = 0
    if dir_ == -1:
        sig = 0
    sigs[ind] = sig

    if dir_ != t_dir:
        ll_max, ll_min = l_max, l_min

    t_dir = dir_

s_df['one_sig'] = sigs

# (s_df['chg'] + 1).cumprod().plot()
(s_df['chg'] * s_df['one_sig'].shift(1) + 1).cumprod().plot(label='stra2')
# plt.show()

# %% 策略2 第二次趋势确认2点时，进行trigger操作。
sig = 0
sigs = np.zeros_like(macd_dir)
s_, t_dir = s[0], macd_dir[0]
ll_max, lll_max, llll_max = max_last[0], max_last[0], max_last[0]
ll_min, lll_min, llll_min = min_last[0], min_last[0], min_last[0]
for ind, (s_, dir_, l_max, l_min) in enumerate(zip(s, macd_dir, max_last, min_last)):

    if dir_ == 1:
        if l_min > lll_min and s_ > ll_max:
            sig = 1
        else:
            sig = 0
    if dir_ == -1:
        sig = 0
    sigs[ind] = sig

    if dir_ != t_dir:
        llll_max, llll_min = lll_max, lll_min
        lll_max, lll_min = ll_max, ll_min
        ll_max, ll_min = l_max, l_min

    t_dir = dir_

s_df['double_sig'] = sigs

(s_df['chg'] * s_df['double_sig'].shift(1) + 1).cumprod().plot(label='stra3')
plt.legend()
plt.show()

# %%
plt.figure(figsize=(16, 8))
plt.title('final result')
(s_df['chg'] + 1).cumprod().plot(label='index')
(s_df['chg'] * sign_long.shift(1) + 1).cumprod().plot(label='strategy')
plt.legend()
plt.show()