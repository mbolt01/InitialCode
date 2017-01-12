# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 16:44:57 2016

@author: mb22
"""
#%%

from lmfit import Model, Parameter
import numpy as np
import matplotlib.pyplot as plt
#%%
def decay(t,n,tau):
    return n*np.exp(-t/tau)
    
#%%

t = np.linspace(0,5,num=1000)
data=decay(t,7,3) + np.random.randn(*t.shape)

#%%

model = Model(decay,independant_vars=['t'])
result = model.fit(data, t=t,
                   n=Parameter(value=10),
                   tau=Parameter(value=1, vary=True))

print(result.values)
print(result.params)

plt.plot(t,data)
plt.plot(t,decay(t=t, **result.values))