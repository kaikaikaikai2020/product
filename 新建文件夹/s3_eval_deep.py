from glob import glob
import pandas as pd
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt

# %% deep eval
deep_list = glob('./deep_result/*.csv')
deep_list.sort()

rets = {}
for path in deep_list:
    date = path[-11:-4]
    df = pd.read_csv(path, index_col=0)
    df.sort_values('pred', ascending=False, inplace=True)
    rets[date] = df.iloc[:35]['real_return'].mean() - df.iloc[-35:]['real_return'].mean()

# %%
(pd.Series(rets).sort_index() + 1 - 0.005).cumprod().plot(figsize=(16, 8))
plt.show()

