# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 15:04:13 2020

@author: adair-9960
"""

pnl = pd.DataFrame()
pnl['Time'] = EG_Test.timestamp
pnl['Trade_Profit'] = 0
pnl['Cost'] = 0
pnl['N_Trade'] = 0

import warnings
warnings.filterwarnings('ignore')
for i in indices:
    if i % 100 == 0:
        print(i)
    EG_Test.process(index=i, transaction_cost=0.001, **action_dict)
    trade_record = EG_Test.record
    if (trade_record is not None) and (len(trade_record) > 0):
        print('value at {}'.format(i))
        trade_record = pd.DataFrame(trade_record)
        trade_cost   = trade_record.groupby('trade_time')['trade_cost'].sum()
        close_cost   = trade_record.groupby('close_time')['close_cost'].sum()
        profit       = trade_record.groupby('close_time')['profit'].sum()
        open_pos     = trade_record.groupby('trade_time')['long_short'].sum()
        close_pos    = trade_record.groupby('close_time')['long_short'].sum() * -1

        pnl['Cost'].loc[pnl['Time'].isin(trade_cost.index)] += trade_cost.values
        pnl['Cost'].loc[pnl['Time'].isin(close_cost.index)] += close_cost.values
        pnl['Trade_Profit'].loc[pnl['Time'].isin(close_cost.index)] += profit.values
        pnl['N_Trade'].loc[pnl['Time'].isin(trade_cost.index)] += open_pos.values
        pnl['N_Trade'].loc[pnl['Time'].isin(close_cost.index)] += close_pos.values

warnings.filterwarnings(action='once')

# Plot the testing result
pnl['PnL'] = (pnl['Trade_Profit'] - pnl['Cost']).cumsum()
plt.plot(pnl['PnL'])
plt.plot(pnl['N_Trade'])
plt.plot(pnl['Time'], pnl['PnL'])

plt.plot(pnl['Time'], pnl['N_Trade'])