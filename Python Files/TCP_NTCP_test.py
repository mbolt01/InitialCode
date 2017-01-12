# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 13:32:54 2016

@author: mb22
"""

## combination of TCP and NTCP calc sharing doses per fraction
## for each created patient specified.

#%%
import importlib
importlib.reload(TCP_NTCP)

#%%

import TCP_NTCP
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use('Qt4Agg')

#%%
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 

#%%
## the TCP function

def TCP_NTCP_calc_plot(d):
    """
    Make the whole thing a function to allow passing of a range of doses.
    All other values remain fixed as interested in the contribution 
    due to dose variaiton (e.g. output)
    """
    TCP_calc = True
    NTCP_calc = True
    
    if TCP_calc == True:
        TCP_results = TCP_NTCP.completeTCPcalc(n=100,
                                      alphabeta_use=3,
                                      alphabeta_sd_use=5,## supply as percent
                                      d=2,
                                      d_shift=d,
                                      d_sd=1,
                                      d_trend=0, # vary the trend value
                                      max_d=100,
                                      dose_of_interest=74,
                                      dose_input = [74],
                                      TCP_input = [0.8],
                                      d_list = None,
                                      #n0 = 200,
                                      weights_input = None)
        plt.plot(TCP_results['nom_doses'],TCP_results['TCP_pop'],
                 color='darkgreen',lw=1,alpha=1,zorder=10)
        plt.plot(TCP_results['dose_input'],TCP_results['TCP_input'],
                 color='green',ls='',marker='o',alpha=0.5)
        for i in range(TCP_results['n']):
            plt.plot(TCP_results['nom_doses'],TCP_results['TCPs'][i],
                     color='green',ls='-',alpha=0.1)
    
    ## should be able to get the list of doses etc from the TCP calc to
    ## then use within the NTCP calc.
    
    if NTCP_calc == True:
        NTCP_results = TCP_NTCP.complete_NTCP_calc(d_data=[45,70, 62, 67, 78, 65],
                                                   ntcp_data=[0.1,0.15,0.1,0.3,0.3, 0.19],
                                                   irrad_perc = 100,
                                                   frac_doses=TCP_results['doses'],
                                                   initial_params_ntcp=None,
                                                   max_dose=100,
                                                   ntcp_params={'td50_1':(80,10),
                                                                'v':(0.5,10),
                                                                'm':(0.51,10),
                                                                'n':(0.3,10)},
                                                   fit_vals = True,
                                                   )
    
        print(NTCP_results['pop_fit'])
        plt.plot(NTCP_results['fit_x'],NTCP_results['fit_y'],
                 color='darkred',lw=1,alpha=1,zorder=10)
        if NTCP_results['d_data'] is not None:
            plt.plot(NTCP_results['d_data'],NTCP_results['ntcp_data'],
                     color='red',ls='',marker='o',alpha=0.5)
        for i in range(TCP_results['n']):
            plt.plot(TCP_results['nom_doses'],NTCP_results['patient_ntcps'][i],
                     color='red',ls='-',alpha=0.1)
    plt.ylim(0,1)
    plt.xlim(0,100)
    plt.title('Population Dose Shift (%): ' + str(d))
    #plt.savefig(('plots/plot_'+str(d)+'.png'),dpi=300,bbox_inches='tight')
    plt.show()
    #plt.clf()
    
    the_doses=TCP_results['doses']
    the_cumdoses=NTCP_results['cum_doses']
    the_ntcps=NTCP_results['patient_ntcps']
    the_tcps = TCP_results['TCPs']
    print("TCP fit: ", TCP_results['tcp_fit'])
    print("NTCP fit: ", NTCP_results['ntcp_fit'])
    return {'the_doses':the_doses,
            'the_cumdoses':the_cumdoses,
            'the_tcps':the_tcps,
            'the_ntcps':the_ntcps,
            }
## for the n patients, the doses are stored in this part fo the dictionary
## 1 patient per row.
#%%

dose_vars = [0,5] # 0 to 20% dose variaiton from the standard.
all_results=[]
for dose in dose_vars:
    print('Dose Variation (%)',dose)
    results = TCP_NTCP_calc_plot(dose)
    #test = all_results.append(TCP_NTCP_calc_plot(dose))
print('Completed All')
