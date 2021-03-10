import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

# %%
data_list = os.listdir('./pred_data')
stock_list = []

# %%
for stk in ['GM', 'F', 'TM', 'TSLA']:
    df = pd.read_csv('./pred_data/att_update_%s.csv' % stk, index_col=0)
    df['ret'] = df['CLOSE'] / df['PRE_CLOSE'] - 1
    df['sign'] = (df['mid'] > df['pred']).astype(int)
    (df['ret'] + 1).cumprod().plot()
    (df['ret'] * df['sign'] + 1).cumprod().plot(label = 'sta')
    plt.title(stk)
    plt.legend()
    plt.show()
