import numpy as np
import pandas as pd

# %%
arr = np.random.random(size=(1000, 16))
df_x = pd.DataFrame(arr, columns=['X_%s' % i for i in range(arr.shape[1])])
df_y = pd.Series(np.random.random(size=(1000,)))


# %% 定义基础的计算符号
def delay_rolling_func(x, delay=10, type='mean'):
    if type == 'mean':
        return pd.Series.rolling(x, delay).mean()
    if type == 'sum':
        return pd.Series.rolling(x, delay).sum()
    if type == 'max':
        return pd.Series.rolling(x, delay).max()
    if type == 'min':
        return pd.Series.rolling(x, delay).min()
    if type == 'std':
        return pd.Series.rolling(x, delay).std()
    if type == 'prod':
        return pd.Series.rolling(x, delay).apply(lambda x: x.prod(), raw=True)
    if type == 'delay':
        return pd.Series.shift(x, delay)
    if type == 'argmax':
        return pd.Series.rolling(x, delay).apply(lambda x: x.argmax(), raw=True)
    if type == 'argmin':
        return pd.Series.rolling(x, delay).apply(lambda x: x.argmin(), raw=True)


def cross_delay_func(x, y, delay=10, type='corr'):
    if type == 'corr':
        return x.rolling(delay).corr(y)
    if type == 'cov':
        return x.rolling(delay).cov(y)


def cross_func(x, y, type='add'):
    if type == 'add':
        return x + y
    if type == 'sub':
        return x - y
    if type == 'mul':
        return x * y
    if type == 'div':
        return x / y


def single_func(x, type='log'):
    if type == 'cube':
        return x ** 3
    if type == 'square':
        return x ** 2
    if type == 'sqrt':
        return np.sign(x) * np.sqrt(x.abs())
    if type == 'log':
        return np.sign(x) * np.log(x.abs())
    if type == 'exp':
        return np.exp((x - x.mean()) / x.std())
    if type == 'sign':
        return np.sign(x)
    if type == 'neg':
        return - x
    if type == 'inv':
        return 1 / x
    if type == 'abs':
        return x.abs()


def signal_func(x):
    pass


# %% generate next 该段用于随机生成一个下一代的特征
def geterate_sign_data(df):
    cols = df.columns
    gen_type = np.random.choice(['one', 'two'])
    delay = np.random.choice(['self', 'delay'])
    if gen_type == 'one' and delay == 'self':
        col = np.random.choice(cols)
        func_type = np.random.choice(['cube', 'square', 'sqrt', 'log',
                                      'exp', 'sign', 'neg', 'inv', 'abs'])
        x_new = single_func(df[col], type=func_type)
        new_type = '%s(%s)' % (func_type, col)

    if gen_type == 'one' and delay == 'delay':
        col = np.random.choice(cols)
        window = np.random.choice([5, 10, 15])
        func_type = np.random.choice(['mean', 'sum', 'std', 'min',
                                      'max', 'prod', 'delay', 'argmax', 'argmin'])
        x_new = delay_rolling_func(df[col], delay=window, type=func_type)
        new_type = '%s(%s,%s)' % (func_type, col, window)

    if gen_type == 'two' and delay == 'self':
        col = np.random.choice(cols, 2)
        func_type = np.random.choice(['add', 'sub', 'mul', 'div'])
        x_new = cross_func(df[col[0]], df[col[1]], type=func_type)
        new_type = '%s(%s,%s)' % (func_type, col[0], col[1])

    if gen_type == 'two' and delay == 'delay':
        col = np.random.choice(cols, 2)
        func_type = np.random.choice(['corr', 'cov'])
        window = np.random.choice([5, 10, 15])
        x_new = cross_delay_func(df[col[0]], df[col[1]], type=func_type, delay=window)
        new_type = '%s(%s,%s,%s)' % (func_type, col[0], col[1], window)

    x_new.name = new_type

    # print(new_type)

    return x_new


# %% 这一段是对上一段的优化。

# 该函数随时生成一个下一代的生成符号
def generate_sign(cols, type='random'):
    gen_list = ['one', 'two']
    delay_list = ['self', 'delay']
    gen_type = np.random.choice(gen_list)
    delay = np.random.choice(delay_list)
    if gen_type == 'one' and delay == 'self':
        col = np.random.choice(cols)
        func_type = np.random.choice(['cube', 'square', 'sqrt', 'log',
                                      'exp', 'sign', 'neg', 'inv', 'abs'])
        gen_type = ('one', 'self', col, func_type)

    if gen_type == 'one' and delay == 'delay':
        col = np.random.choice(cols)
        window = np.random.choice([5, 10, 15])
        func_type = np.random.choice(['mean', 'sum', 'std', 'min',
                                      'max', 'prod', 'delay', 'argmax', 'argmin'])
        gen_type = ('one', 'delay', col, window, func_type)

    if gen_type == 'two' and delay == 'self':
        col = np.random.choice(cols, 2).tolist()
        func_type = np.random.choice(['add', 'sub', 'mul', 'div'])
        gen_type = ('two', 'self', col, func_type)

    if gen_type == 'two' and delay == 'delay':
        col = np.random.choice(cols, 2).tolist()
        func_type = np.random.choice(['corr', 'cov'])
        window = np.random.choice([5, 10, 15])
        gen_type = ('two', 'delay', col, window, func_type)

    return gen_type

# 这一段函数接受下一代产生的符号，进行运算。并标明计算方法
def generate_data_by_sign(signs, df):
    gen_type, delay = signs[0], signs[1]

    if gen_type == 'one' and delay == 'self':
        col, func_type = signs[2], signs[3]
        x_new = single_func(df[col], type=func_type)
        new_type = '%s(%s)' % (func_type, col)

    if gen_type == 'one' and delay == 'delay':
        col, window, func_type = signs[2], signs[3], signs[4]
        x_new = delay_rolling_func(df[col], delay=window, type=func_type)
        new_type = '%s(%s,%s)' % (func_type, col, window)

    if gen_type == 'two' and delay == 'self':
        col, func_type = signs[2], signs[3]
        x_new = cross_func(df[col[0]], df[col[1]], type=func_type)
        new_type = '%s(%s,%s)' % (func_type, col[0], col[1])

    if gen_type == 'two' and delay == 'delay':
        col, window, func_type = signs[2], signs[3], signs[4]
        x_new = cross_delay_func(df[col[0]], df[col[1]], type=func_type, delay=window)
        new_type = '%s(%s,%s,%s)' % (func_type, col[0], col[1], window)

    x_new.name = new_type
    return x_new

# 该段函数用来生成不重复的随机产生的下一代计符号
def generate_norepeat_signs(cols, n=100):
    signs_list = []
    while 1:
        signs = generate_sign(cols)
        if signs not in signs_list:
            signs_list.append(signs)
        if len(signs_list) == n:
            break
    return signs_list

# 该段函数用来验证下一代的特征的表现
def correlation_filter(df_x, df_y, chs_number=100):
    # 计算每个特征的ic表现
    ic = df_x.corrwith(df_y, method='spearman').drop_duplicates().dropna()
    sort_cols = ic.abs().sort_values(ascending=False)
    df_sort = df_x[sort_cols.index]
    corr_df = df_sort.corr()
    corr_df.replace(-1, 1, inplace=True)
    cols_continue = sort_cols.index.tolist()

    # 进行相似特征的删除，确保生成的子代
    # 1. 特征与特征之间没有强的相关。 2. 父代与子代容易出现较大相关性，可以进行类似家族竞争的结果。
    flag = 0
    for col in cols_continue:
        flag += 1
        sub = corr_df[(corr_df[col].abs() > 0.7) & (corr_df[col] != 1)]
        if len(sub.index) != 0:
            for sub_col in sub.index:
                cols_continue.remove(sub_col)
        corr_df = corr_df.loc[cols_continue, cols_continue]
        if flag > chs_number:
            break
    return corr_df.index.tolist()[:chs_number]


# %%
# sign_list = generate_norepeat_signs(df_x.columns, 200)
#
# signs = generate_sign(df_x.columns)
# # print(signs)
# x_new = generate_data_by_sign(df=df_x, signs=signs)

# %% 这段函数是多线程多进程来提高计算效率的代码
# if __name__ == '__main__':
#     from multiprocessing import Pool
#     from functools import partial
#
#     # p = Pool(8)
#     for i in range(4):
#         sign_list = generate_norepeat_signs(df_x.columns, 300)
#         # params_list = [(df_x, signs) for signs in sign_list]
#         with Pool(8) as p:
#             sub_xs = p.map(partial(generate_data_by_sign, df = df_x), sign_list)
#
#         sub_dfx = pd.concat(sub_xs, axis=1)
#         sub_dfx = sub_dfx.loc[:, ~sub_dfx.columns.duplicated()]
#
#         df_x = pd.concat([df_x, sub_dfx], axis=1)
#         df_x = df_x.loc[:, ~df_x.columns.duplicated()]
#
#         # avoid nan value
#         df_x = df_x[df_x.columns[df_x.isna().mean() < 0.2]]
#
#         ic = df_x.corrwith(df_y, method='spearman').drop_duplicates().dropna()
#         # get ic abs top
#         chs_cols = ic.abs().sort_values(ascending=False).iloc[:100]
#         print(chs_cols)
#         df_x = df_x[chs_cols.index]

# %%
if __name__ == '__main__':
    sign_list = generate_norepeat_signs(df_x.columns, 300) # 首先生成不重复的300个子代生成的符号
    print(sign_list[:10])
    # 通过符号来进行计算
    sub_xs = [generate_data_by_sign(df=df_x, signs=signs) for signs in sign_list]
    # 将所有生成的子代整合成一个大的数据集
    sub_dfx = pd.concat(sub_xs, axis=1)
    sub_dfx = sub_dfx.loc[:, ~sub_dfx.columns.duplicated()]

    df_x = pd.concat([df_x, sub_dfx], axis=1)
    df_x = df_x.loc[:, ~df_x.columns.duplicated()]

    print(df_x.head())

    # 将所有特征送入进行筛选，选择出最优的特征作为下一次迭代的父代。
    chs_cols = correlation_filter(df_x, df_y)

    # ic = df_x.corrwith(df_y).drop_duplicates().dropna()
    # # get ic abs top
    # chs_cols = ic.abs().sort_values(ascending=False).iloc[:100]
    print(chs_cols)
    df_x = df_x[chs_cols]

# %%
# df_scale = (df_x - df_x.mean()) / df_x.std()

# %%


# %% filter the corr
