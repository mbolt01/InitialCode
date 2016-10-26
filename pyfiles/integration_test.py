# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 15:12:29 2016

@author: mb22
"""

import scipy as sp
import numpy as np

#%%

### option 1 = supply the integral to solve
def integrand(x):
    return np.exp(-0.5*x**2)

def int_test(a,b):
    # a and b are the lower and upper bounds of the integration
    ## this shoudl work if limits are sensible (u < 36 works fine)
    ## this should be the case for me otherwise the result will be 100% NTCP
    return sp.integrate.quad(integrand,a,b,limit=500)
    
### option 2 = integral is already solved and supply the limits to this
def alt_calc(x):
    """
    This is the solution to the required integral.
    This seems more reliable as I do not know the required limits.
    """
    return sp.sqrt(sp.pi)*sp.special.erf(x/sp.sqrt(2))/sp.sqrt(2)
    
def alt_int_test(a,b):
    return alt_calc(b)-alt_calc(a)
    
### Show both give the same results    
def both(a,b):
    return (int_test(a,b)[0],alt_int_test(a,b))
    
#%%

import matplotlib.pyplot as plt
allres = []
alli=[]

for i in range(500):
    result = both(-np.inf,i)
    dif = result[0]-result[1]
    alli.append(i)
    allres.append(dif)

plt.plot(alli,allres)
#plt.ylim(-0.0000000001,0.0000000001)
## randomly fails at i=21, and any value greater than 35
plt.show()