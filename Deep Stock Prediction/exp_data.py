import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# %%
df = pd.read_csv('./stk_data/TSLA.O.csv', index_col=0)

df['CLOSE'].plot()
plt.show()

# %%