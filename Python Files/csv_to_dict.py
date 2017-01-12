# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 12:39:01 2016

@author: mb22
"""
#%%

import csv

FILE = csv.reader(open(r'C:\Users\mb22\OneDrive\PhD\Work\results2.csv'))

#%%

lines = [row for row in FILE]
params = lines[0]
vals = lines[1]

results_dict = dict(zip(params,vals))

#%%

params = next(FILE)
vals = next(FILE)

results_dict = dict(zip(params,vals))

#%%

cnr_result = results_dict['CNR Air_CNR']
recon_diam = float(results_dict['RECONSTRUCTION_DIAMETER'])
print(cnr_result)
print(recon_diam)
print(type(cnr_result))
print(type(recon_diam))


#%%
