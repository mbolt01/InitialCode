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


#CHHIP: node negative T1b-T3a localised PCa
# with risk of seminal vesical involvement ≤30%

all_test = TCP_full(k=5,
                    TCP_input=88.3,
                    repeats=5,
                    n=500,
                    alphabeta_use=3,
                    alphabeta_sd_use=0.5,
                    d=2,
                    d_shift=0,
                    d_sd=0.5,
                    d_trend=0,
                    max_d=100,
                    dose_of_interest=74,
                    save_name="Results16June16-CHHIP_74-Trend=0-ab_sdvar")

print("Completed....")

#%%
## load CSV file into pandas for analysis
          
TCP_data = pd.read_csv(os.getcwd()+"\\"+all_test[0]) # read in file saved after TCP simulation.

print(all_test[0])
TCP_data.describe()

#%%
