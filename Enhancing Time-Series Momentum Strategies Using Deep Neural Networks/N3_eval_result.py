import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# %%
n_stk = 3
for alg_type in ['binary', 'reg', 'avg', 'shp']:
    eval_df = pd.read_csv('./eval_data/%s.csv'%alg_type)
    eval_df.set_index(['tradeDate', 'code'], inplace=True)
    df_sig = eval_df['pred_val'].unstack()
    df_ret = eval_df['eval_ret'].unstack()
    ret_dict = {}
    trade_date = eval_df.index.get_level_values('tradeDate').unique()
    for i in tqdm(trade_date):
        sub_df = eval_df.loc[i]
        chs_stk = sub_df.abs().sort_values('pred_val', ascending=False).index[:n_stk]
        ret_dict[i] = (sub_df.loc[chs_stk, 'eval_ret'] *
                       sub_df.loc[chs_stk, 'pred_val'].apply(lambda x: np.sign(x))).mean()


    (pd.Series(ret_dict) + 1).cumprod().plot(label = alg_type)

(df_ret.mean(axis = 1) + 1).cumprod().plot(label = 'hold long')
plt.legend()
plt.show()
