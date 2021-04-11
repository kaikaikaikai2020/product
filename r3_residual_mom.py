import pandas as pd
from glob import glob
from tqdm import tqdm
import numpy as np
import statsmodels.api as sm

# %%
stk_list = glob('./fac_data/*.csv')


# %%
def get_mom_data(path, window_size=36, mom_step=12):
    df = pd.read_csv(path, index_col=0)
    df['mom_residual'] = np.nan
    df['mom_convent'] = np.nan

    month_list = df.index.tolist()
    n = df.shape[0]
    for i in range(n - window_size):
        sub_df = df.iloc[i: i + window_size]
        sub_x = sub_df[['CHGPct', 'smb', 'hml']]
        sub_y = sub_df['chg_pct']
        model = sm.OLS(sub_y, sm.add_constant(sub_x)).fit()
        df.loc[month_list[i + window_size], 'mom_residual'] = (model.resid.iloc[-mom_step:-2] + 1).prod() - 1
        df.loc[month_list[i + window_size], 'mom_convent'] = (sub_y.iloc[-mom_step:-2] + 1).prod() - 1
    return df.dropna()

# %%
import multiprocessing

if __name__ == '__main__':
    s_list = []
    num_core = 8
    file_step = 16
    n = int(len(stk_list) / 16 + 1)

    for i in tqdm(range(n)):
        with multiprocessing.Pool(num_core) as pool:
            s = pool.map(get_mom_data, stk_list[i * file_step: (i + 1) * file_step])
            s_list.extend(s)

    df = pd.concat(s_list)
    df.to_csv('./moment.csv')

# %%
