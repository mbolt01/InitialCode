# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 10:21:07 2016

@author: mb22
"""

import TCP
import matplotlib.pyplot as plt
import importlib
importlib.reload(TCP)

TCP = 88
n = 500
alphabeta_use = 3
alphabeta_sd_use = 0.5
d = 2
d_shift = 0
d_sd = 0.5
d_trend = 0
max_d = 100
dose_of_interest = 74
TCP_input = 88
   


xyz =  TCP.calc_dif_sq(x,TCP, n, alphabeta_use, alphabeta_sd_use,d,d_shift,d_sd,d_trend,max_d,dose_of_interest,TCP_input)
    