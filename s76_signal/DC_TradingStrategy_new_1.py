# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 14:00:53 2020

@author: Asus
"""




import matplotlib.pyplot as plt
from yq_toolsS45 import get_MktEqudAdjAfGet_update 
df = get_MktEqudAdjAfGet_update('000002','2010-01-01','2020-12-18',
                      'tradeDate as Date,openPrice as Open, highestPrice as High,lowestPrice as Low,closePrice as Close, turnoverVol as Volume')
#df.set_index(['Date'],drop=True,inplace=True)
df.sort_values(by=['Date'])
#print(df)


index_list=df['Close']
#--------初始化变量---------------
upturn_event = True
p_h = p_l = index_list[0]
dc_range = []
os_range = []
dc_direction = []
r=0
thresh=0.01
sum1=100000
osup=0
osdown=0

#---------------寻找DCC和OS点----------
for i in range (len(index_list)):
    p_t = index_list[i]
    if upturn_event:
        if p_t <= p_h * (1 - thresh): #最高点下跌一个theta
            upturn_event = False  #更改up标志位
            p_l = p_t #低点为当前点

            dc_range.append(index_list[i])
            dc_direction.append('S')
            os_range.append(index_list[osup])
            #table_wt.write(i,8,index_list[i])
            #excel.save("SH2.xlsx")
            #table_wt.write(i,10,index_list[osup])
            #excel.save("SH2.xlsx")
            
        else:
            
            if p_h < p_t:
                p_h = p_t
                osup=i

    else: # if we're in a downturn
        if p_t >= p_l * (1 + thresh):
            upturn_event = True
            p_h = p_t
            dc_range.append(index_list[i])
            dc_direction.append('B')
            os_range.append(index_list[osdown])
            #table_wt.write(i,9,index_list[i])
            #excel.save("SH2.xlsx")
            #table_wt.write(i,11,index_list[osdown])
            #excel.save("SH2.xlsx")
        
        else:
            if p_l > p_t:
                p_l = p_t
                osdown=i


print (dc_range)
print(dc_direction)
print(os_range)
temp = []
#----------profit----------------
for x in range (1,len (dc_range)-1,1):
    
    if x%2==1:
        profit=(dc_range[x+1]-dc_range[x])/dc_range[x]
        sum1=sum1*(1+profit)
    else:
        profit=(dc_range[x]-dc_range[x+1])/dc_range[x]
        sum1=sum1*(1+profit)
    temp.append(sum1)
print (sum1)
plt.plot(temp)

