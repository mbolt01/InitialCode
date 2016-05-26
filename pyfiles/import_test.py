# -*- coding: utf-8 -*-
"""
Created on Thu May 26 16:12:28 2016

@author: mb22
"""

#%%

from TCP import *

#%%

## perform the analysis and save as CSV file.

print("Started....")

all_test = TCP_full(k=5,
                    TCP_input=75,
                    repeats=6,
                    n=80,
                    alphabeta_use=3,
                    alphabeta_sd_use=0.2,
                    d=2,
                    d_shift=-1,
                    d_sd=0.6,
                    d_trend=0.5,
                    max_d=100,
                    dose_of_interest=74,
                    save_name="Results1-Trend=0")

print("Completed....")

#%%
## load CSV file into pandas for analysis                  
TCP_data = pd.read_csv(os.getcwd()+"\\"+all_test[0]) # read in file saved after TCP simulation.

print(all_test[0])
TCP_data.describe()

#%%
