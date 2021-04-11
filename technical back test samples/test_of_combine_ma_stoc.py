# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 19:09:28 2021

@author: ASUS
"""
margin = 0.001
def signal(Data, what, moving_average, indicator, buy, sell):
    
    for i in range(len(Data)):
        if Data[i, what] > Data[i, moving_average] and Data[i, what] -  Data[i, moving_average] <= margin and Data[i - 1, what] -  Data[i - 1, moving_average] > margin and Data[i, indicator] <= lower_barrier and Data[i - 1, indicator] > lower_barrier:
            Data[i, buy] = 1
    
        if Data[i, what] < Data[i, moving_average] and Data[i, moving_average] -  Data[i, what] <= margin and Data[i - 1, moving_average] -  Data[i - 1, what] > margin and Data[i, indicator] >= upper_barrier and Data[i - 1, indicator] < upper_barrier:
            Data[i, sell] = -1