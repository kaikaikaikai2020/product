import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# %%
class minmax_process:

    def __init__(self, x, macd_x):
        self.macd_dir = macd_x
        self.s_df = x[['closeIndex']].copy()
        self.s = x['closeIndex'].copy()

    def spilt_dir(self):
        self.s_df['direction'] = self.macd_dir
        pre_direction = self.s_df['direction'].shift(1)
        flag = 0
        flags = []
        for i, j in zip(self.macd_dir, pre_direction):
            if i != j:
                flag += 1
            flags.append(flag * i)
        self.s_df['flag_part'] = flags
        max_dict = self.s_df.groupby(flags)['closeIndex'].max().to_dict()
        min_dict = self.s_df.groupby(flags)['closeIndex'].min().to_dict()
        self.s_df['dir_max'] = [max_dict[i] for i in flags]
        self.s_df['dir_min'] = [min_dict[i] for i in flags]

    def generate_dir_df(self):
        dir_df = self.s_df.groupby('flag_part')[['dir_min', 'dir_max']].mean()
        dir_df['dir_flag'] = dir_df.index
        dir_df = dir_df.iloc[np.argsort(dir_df.index.map(np.abs)), :]
        self.dir_df = dir_df

    def update_dir(self):
        min_ser, max_ser, dir_ser = self.dir_df['dir_min'].values, \
                                    self.dir_df['dir_max'].values, \
                                    self.dir_df['dir_flag'].values
        new_dir = np.zeros_like(dir_ser)
        flag = True

        for ind, (min_, max_, dir_) in enumerate(zip(min_ser, max_ser, dir_ser)):
            chg_flag = False
            if flag:
                t_min, t_max, t_dir = min_, max_, dir_
                flag = False
            if dir_ > 0:
                if min_ < t_min and max_ < t_max:
                    new_ = -dir_
                    chg_flag = True
            if dir_ < 0:
                if max_ > t_max and min_ > t_min:
                    new_ = -dir_
                    chg_flag = True

            # 记录临时值
            t_min, t_max, t_dir = min_, max_, dir_

            if chg_flag:
                new_dir[ind] = np.sign(new_)
            else:
                new_dir[ind] = np.sign(dir_)

        self.dir_df['update_dir'] = new_dir

    def update_except(self):
        flag = True
        t_exp = 1
        min_ser, max_ser, new_dir = self.dir_df['dir_min'].values, \
                                    self.dir_df['dir_max'].values, \
                                    self.dir_df['update_dir'].values
        last_max, last_min = np.zeros_like(new_dir), np.zeros_like(new_dir)
        exp_ser = np.zeros_like(new_dir)
        for ind, (min_, max_, dir_) in enumerate(zip(min_ser, max_ser, new_dir)):
            chg_flag = False
            if flag:
                t_min, t_max, t_dir = min_, max_, dir_
                flag = False

            if t_exp == -1:
                cd1 = t_dir * dir_ == -1
                cd2 = t_dir == 1 and t_max < max_
                cd3 = t_dir == -1 and t_min > min_
                chg_flag = True
                if cd1 or cd2 or cd3:
                    exp = 1
                else:
                    exp = -1
                # print('1', exp)
            elif t_dir == dir_:
                cd1 = dir_ == 1 and t_min > min_
                cd2 = dir_ == -1 and t_max < max_
                chg_flag = True
                if cd1 or cd2:
                    exp = -1
                else:
                    exp = 1
                # print('2', exp)

            if t_dir == dir_:
                t_min, t_max, t_dir = np.min([min_, t_min]), np.max([max_, t_max]), dir_
            else:
                t_min, t_max, t_dir = min_, max_, dir_

            if chg_flag:
                exp_ser[ind] = exp
                t_exp = exp
            else:
                exp_ser[ind] = 1

            last_min[ind] = t_min
            last_max[ind] = t_max

            t_dir = dir_

        self.dir_df['exp'] = exp_ser
        # update max min ser
        self.dir_df['dir_max'] = last_max
        self.dir_df['dir_min'] = last_min
        self.dir_df['stat'] = self.dir_df['dir_flag'].apply(np.sign) * self.dir_df['exp']

    def update_minmax(self):
        flag = True
        min_ser, max_ser, stat_ser = self.dir_df['dir_min'].values, \
                                     self.dir_df['dir_max'].values, \
                                     self.dir_df['stat'].values
        last_max, last_min = np.zeros_like(stat_ser), np.zeros_like(stat_ser)
        for ind, (min_, max_, stat_) in enumerate(zip(min_ser, max_ser, stat_ser)):
            if flag:
                t_min, t_max, t_stat = min_, max_, stat_
                l_min, l_max = t_min, t_max
                flag = False

            if t_stat == 1:
                if t_max < max_:
                    t_max = max_

                if stat_ == -1:
                    l_max = t_max
                    t_min = np.min([t_max, l_min, min_])

            if t_stat == -1:
                if t_min > min_:
                    t_min = min_

                if stat_ == 1:
                    l_min = t_min
                    t_max = np.max([l_min, t_max, max_])

            last_max[ind] = l_max
            last_min[ind] = l_min

            t_min, t_max, t_stat = min_, max_, stat_
            l_min, l_max = t_min, t_max

        self.dir_df['last_max'] = last_max
        self.dir_df['last_min'] = last_min

    def fit(self):
        self.spilt_dir()
        self.generate_dir_df()
        self.update_dir()
        self.update_except()
        self.update_minmax()
        self.s_df['stat'] = self.s_df['flag_part'].map(self.dir_df['stat'].to_dict())
        self.s_df['last_max'] = self.s_df['flag_part'].map(self.dir_df['last_max'].to_dict())
        self.s_df['last_min'] = self.s_df['flag_part'].map(self.dir_df['last_min'].to_dict())
        self.s_minmax = self.s_df.apply(lambda x: x['dir_max'] if x['stat'] == 1 else x['dir_min'], axis=1)

# %%
if __name__ == '__main__':
    from alg_macd_direction import macd_direction
    df = pd.read_csv('./data/hs_300.csv', index_col=0)
    s = df['closeIndex']
    macd_dir = macd_direction().integal_macd(s)


    result = minmax_process(df, macd_dir)
    result.fit()
    plt.figure(figsize=(16, 5))
    plt.plot(result.s[1500:2000])
    plt.plot(result.s_minmax[1500:2000])
    plt.show()
# %%
# s_df = df[['closeIndex']].copy()
# s_df['direction'] = macd_dir
# pre_direction = s_df['direction'].shift(1)
#
# # %%
# flag = 0
# flags = []
# for i, j in zip(s_df['direction'], pre_direction):
#     if i != j:
#         flag += 1
#     flags.append(flag * i)
# s_df['flag_part'] = flags
# max_dict = s_df.groupby(flags)['closeIndex'].max().to_dict()
# min_dict = s_df.groupby(flags)['closeIndex'].min().to_dict()
# s_df['dir_max'] = [max_dict[i] for i in flags]
# s_df['dir_min'] = [min_dict[i] for i in flags]
#
# # %% 更新dir
# dir_df = s_df.groupby('flag_part')[['dir_min', 'dir_max']].mean()
# dir_df['dir_flag'] = dir_df.index
# dir_df = dir_df.iloc[np.argsort(dir_df.index.map(np.abs)), :]
#
# # %%
# min_ser, max_ser, dir_ser = dir_df['dir_min'].values, dir_df['dir_max'].values, dir_df['dir_flag'].values
# new_dir = np.zeros_like(dir_ser)
# flag = True
#
# for ind, (min_, max_, dir_) in enumerate(zip(min_ser, max_ser, dir_ser)):
#     chg_flag = False
#     if flag:
#         t_min, t_max, t_dir = min_, max_, dir_
#         flag = False
#     if dir_ > 0:
#         if min_ < t_min and max_ < t_max:
#             new_ = -dir_
#             chg_flag = True
#     if dir_ < 0:
#         if max_ < t_max and min_ > t_min:
#             new_ = -dir_
#             chg_flag = True
#
#     # 记录临时值
#     t_min, t_max, t_dir = min_, max_, dir_
#
#     if chg_flag:
#         new_dir[ind] = np.sign(new_)
#     else:
#         new_dir[ind] = np.sign(dir_)
#
# # %%
# dir_df['update_dir'] = new_dir
#
# # %% 更新exp参数,stat
# flag = True
# t_exp = 1
# exp_ser = np.zeros_like(new_dir)
# for ind, (min_, max_, dir_) in enumerate(zip(min_ser, max_ser, new_dir)):
#     chg_flag = False
#     if flag:
#         t_min, t_max, t_dir = min_, max_, dir_
#         flag = False
#
#     if t_exp == -1:
#         cd1 = t_dir * dir_ == -1
#         cd2 = t_dir == 1 and t_max < max_
#         cd3 = t_dir == -1 and t_min > min_
#         chg_flag = True
#         if cd1 or cd2 or cd3:
#             exp = 1
#         else:
#             exp = -1
#         # print('1', exp)
#     elif t_dir == dir_:
#         cd1 = dir_ == 1 and t_min > min_
#         cd2 = dir_ == -1 and t_max < max_
#         chg_flag = True
#         if cd1 or cd2:
#             exp = -1
#         else:
#             exp = 1
#         # print('2', exp)
#
#     if t_dir == dir_:
#         t_min, t_max, t_dir = np.min([min_, t_min]), np.max([max_, t_max]), dir_
#     else:
#         t_min, t_max, t_dir = min_, max_, dir_
#
#     if chg_flag:
#         exp_ser[ind] = exp
#         t_exp = exp
#     else:
#         exp_ser[ind] = 1
#
# # %%
# dir_df['exp'] = exp_ser
# dir_df['stat'] = dir_df['dir_flag'].apply(np.sign) * dir_df['exp']
#
# # %% minmax processing
# flag = True
# stat_ser = dir_df['stat'].values
# last_max, last_min = np.zeros_like(stat_ser), np.zeros_like(stat_ser)
# for ind, (min_, max_, stat_) in enumerate(zip(min_ser, max_ser, stat_ser)):
#     if flag:
#         t_min, t_max, t_stat = min_, max_, stat_
#         l_min, l_max = t_min, t_max
#         flag = False
#
#     if t_stat == 1:
#         if t_max < max_:
#             t_max = max_
#
#         if stat_ == -1:
#             l_max = t_max
#             t_min = np.min([t_max, l_min, min_])
#
#     if t_stat == -1:
#         if t_min > min_:
#             t_min = min_
#
#         if stat_ == 1:
#             l_min = t_min
#             t_max = np.max([l_min, t_max, max_])
#
#     last_max[ind] = l_max
#     last_min[ind] = l_min
#
#     t_min, t_max, t_stat = min_, max_, stat_
#     l_min, l_max = t_min, t_max
#
# # %%
# s_df['stat'] = s_df['flag_part'].map(dir_df['stat'].to_dict())
# s_minmax = s_df.apply(lambda x: x['dir_max'] if x['stat'] == 1 else x['dir_min'], axis=1)
#
# # %%
# plt.plot(s[1000:1500])
# plt.plot(s_minmax[1000:1500])
# plt.show()
