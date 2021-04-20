import pandas as pd
import numpy as np
import os
from glob import glob
import zipfile
from tqdm import tqdm

# %%
ms = glob('./data/*/*.zip') + glob('./data/*/*/*.zip')
ms = [i for i in ms if 'fenbi' not in i]
if not os.path.exists('./unzip_data/'): os.mkdir('./unzip_data/')


# %% unzip data
def zip_file_func(path, to_path):
    zip_file = zipfile.ZipFile(path)
    a_name = zip_file.namelist()
    for i in a_name:
        if '1min' in i:
            break
    zip_file.extract(i, to_path)


#for path in tqdm(ms):
#    try:
#        zip_file_func(path, './unzip_data/')
#    except:
#        print(path)

# %% processing data
cols = ['code', 'time', 'open', 'high', 'low', 'close', 'vol', 'value', 'amount']
ms = glob('./unzip_data/*.csv') + glob('./unzip_data/*/*.csv')


# %%
def processing_data(path):
    try:
        df = pd.read_csv(path, encoding='gbk', header=1)
        df.columns = cols
        date = df['time'].loc[0].split(' ')[0]
        df_day = pd.pivot_table(df, index=['code'], values=['open', 'high', 'low', 'close'],
                                aggfunc={'open': lambda x: x.iloc[0],
                                         'high': lambda x: x.max(),
                                         'low': lambda x: x.min(),
                                         'close': lambda x: x.iloc[-1]})

        df_twap = df_day.groupby('code').apply(lambda x: x[['open', 'high', 'low', 'close']].mean().mean())
        df_day['twap'] = df_twap
        df_day['arpp'] = (df_day['twap'] - df_day['low']) / (df_day['high'] - df_day['low'])
        return df_day, date
    except:
        print(date)
        return None, None

# %%
#for path in tqdm(ms):
#    df, date = processing_data(path)
#    df.to_csv('./date_data/%s.csv' % date)
#    print(date)


# %%
if __name__ == '__main__':
    from multiprocessing import Pool

    if not os.path.exists('./date_data/'): os.mkdir('./date_data/')

    n_step = 16
    n_round = int(len(ms) / n_step) + 1
    print(n_round)
    for i in tqdm(range(n_round)):
        with Pool(8) as pool:
            s = pool.map(processing_data, ms[i * n_step: (i + 1) * n_step])
            for (df, date) in s:
                if date:
                    df.to_csv('./date_data/%s.csv' % date)
