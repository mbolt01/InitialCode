# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 09:31:24 2017

@author: mb22
"""

## Run TCP/NTCP simulation for supplied parameters.
## save the output as a pickel for later access

import TCP_NTCP
#%%
import matplotlib.pyplot as plt
#%%
import numpy as np
#%%
import pickle
#%%
import sys
import random

import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 

#%%
import importlib
importlib.reload(TCP_NTCP)

#%%

## get results from the TCP simulation

save = False

my_sds = [0, 0.5,
          1, 1.5,
          2, 2.5,
          3, 3.5,
          4, 4.5,
          5, 5.5,
          6, 6.5,
          7, 7.5,
          8, 8.5,
          9, 9.5,
          10] ## demonstrate the impact of the daily output variation.
          
my_sds = [5]


my_outputs = [-2,-1,0,1,2] ## use this to model single patient on multiple linacs
#my_outputs=[-2]

for my_sd in my_sds:
    TCP_results = TCP_NTCP.completeTCPcalc(n=100,
                                          alphabeta_use=10,
                                          alphabeta_sd_use=0,
                                          d=65/30,
                                          d_shift=0, # no shift as start each treatment as if perfect
                                          d_sd=my_sd,
                                          d_trend=0, # vary the trend value
                                          max_d=100,
                                          dose_of_interest=65,
                                          dose_input = [65],
                                          TCP_input = [0.66],
                                          d_list = None,
                                          n0 = None,
                                          weights_input = None)
    
    ##%%
    
    ## get results from the NTCP simulation
    
    ## can specify different NTCP doses if required.
    ## usually want the same as modelling same patient(s)
    #ntcp_doses = TCP_results['doses']*2
    
    NTCP_results = TCP_NTCP.complete_NTCP_calc(d_data=[35,59],
                                               ntcp_data=[0.29,0.83],
                                               irrad_perc = 100, ## scaling factor?
                                               frac_doses=TCP_results['doses'],
                                               #frac_doses=ntcp_doses, ## can specify different set of doses
                                               #initial_params_ntcp=None,
                                               max_dose=TCP_results['max_d'],
                                               #ntcp_params={'td50_1':(58.2,1.92),
                                               #             'v': None,#(0.08,10),
                                               #             'm':(0.28,37.3),
                                               #             'n':(0.14,16.43)}
                                               initial_params_ntcp = [31.4,
                                                                      0.7,
                                                                      0.3,
                                                                      1],
                                               ntcp_params={'td50_1':(31.4,0),
                                                            'v': (0.7,0),#(0.08,10),
                                                            'm':(0.35,0),
                                                            'n':(1,0)}, #(1,0)
                                               fit_vals=False)
    
    ##%%
    
    ## extract the data from all patients at a set dose point
    ## i.e. get a column from the numpy results array to then do stata

    my_dose_tcp = 65 ## choose dose of interest (usually nominal dose for tcp)
    my_dose_ntcp = 30 ## different dose of interest for ntcp
    
    ## get the position in the array of nominal doses for use.
    index_tcp, value_tcp = TCP_NTCP.closest_val(TCP_results['nom_doses'],my_dose_tcp)
    index_ntcp, value_ntcp = TCP_NTCP.closest_val(TCP_results['nom_doses'],my_dose_ntcp)
    #print('Dose Value:', value, 'Gy, Array Index:',index)
    
    ## get the column as a variable so can then do stats on these only.
    the_TCP= TCP_results['TCPs'][:,index_tcp]
    the_NTCP = NTCP_results['patient_ntcps'][:,index_ntcp]
    
    ## get some stats and store in dict using function
    def stat_dict(results):
        """ Calculate and store statistical results at dose of interest.
        The list of values of TCP/NTCP at the dose of interest should be supplied
        """
        
        stats = {'mean':np.mean(results),
                 'sd':np.std(results),
                 'perc5':np.percentile(results,5),
                 'perc95':np.percentile(results,95),
                 'median':np.median(results),
                 'max':np.max(results),
                 'min':np.min(results),
                 'perc25':np.percentile(results,25),
                 'perc75':np.percentile(results,75),
                 
                 'sdperc':100*np.std(results)/np.mean(results),
                 'range95_5':np.percentile(results,95)-np.percentile(results,5),
                 'range95_5_perc':100*(np.percentile(results,95)-np.percentile(results,5))/np.mean(results),
                 }
        
        return stats
        
    tcp_stats = stat_dict(the_TCP)
    ntcp_stats = stat_dict(the_NTCP)
    
    
    ##%%
    
    ## save the results in a pickle file
    
    #all_results = {'TCP':TCP_results,'NTCP':NTCP_results}
    
    the_name = 'results/headneck/headneck_'+str(my_sd)+'perc.pkl'
    
    if save == True:
        with open(the_name, 'wb') as f:
            ## note the results are stored in a dict for ease of later access
            ## do this rather than adding each in turn as then order is important.
            pickle.dump({'TCP':TCP_results,'NTCP':NTCP_results,
                         'TCPstats':tcp_stats,
                         'NTCPstats':ntcp_stats}, f)
        
        print('Completed:', my_sd)
    else:
        print('Completed:', my_sd, ' (Result not saved)')
    
        
#%%
## plot the results of the simulation

TCP_NTCP.plot_TCP_NTCP(resultsTCP=TCP_results, resultsNTCP=NTCP_results,
              n=0, colors={'TCP':'green','NTCP':'red'},
              pop_plot=True, plot_percentiles=[5,95], show_percentiles=True,
              TCP=False, show_legend=False)

## plot vertical line at my_dose point where stats are relevant
#plt.axvline(my_dose, color='black', alpha=0.8, lw=0.5)

## plot line to indicate the points of interest to show 95th percentle range
#plt.plot([my_dose,my_dose],[np.percentile(the_NTCP,5), np.percentile(the_NTCP,95)],
#         marker='', color='orange', zorder=10, lw=5, alpha=0.8)
#plt.plot([my_dose,my_dose],[np.percentile(the_TCP,5), np.percentile(the_TCP,95)],
#         marker='', color='orange', zorder=10, lw=5, alpha=0.8)

    
#%%
## load the results back from the saved files to then plot/analyse

meansTCP = {}
meansNTCP = {}

sdsTCP = {}
sdsNTCP = {}

range95TCP = {}
range95NTCP = {}

for my_sd in my_sds:
    the_name = 'results/headneck/headneck_'+str(my_sd)+'perc.pkl'
    
    with open(the_name, 'rb') as f:
        ## load the file then extract the different parts of the dict
        ## loading from dict removes the need to know the order.
        loaded = pickle.load(f)
        TCP=loaded['TCP']
        NTCP=loaded['NTCP']
        TCPstats=loaded['TCPstats']
        NTCPstats=loaded['NTCPstats']
        
    meansTCP[my_sd] = TCPstats['mean']
    meansNTCP[my_sd] = NTCPstats['mean']
    
    sdsTCP[my_sd] = TCPstats['sdperc']
    sdsNTCP[my_sd] = NTCPstats['sdperc']

    range95TCP[my_sd] = TCPstats['range95_5_perc']
    range95NTCP[my_sd] = NTCPstats['range95_5_perc']

#%%
## plot some of the results for range of dose variations
to_plots = [#('TCP mean',meansTCP),
            #('NTCP mean',meansNTCP),
            ('Range TCP (95th Percentiles)',range95TCP),
            ('Range NTCP (95th Percentiles)',range95NTCP)]
            
#print(to_plots[0][0])
for i in range(len(to_plots)):
    data_label = to_plots[i][0]
    data = to_plots[i][1]

    x, y = zip(*data.items())
    plt.plot(x,y, 'o', ms=4, label=data_label)
plt.ylabel('Value (%)')
plt.xlabel('Dose Standard Deviation (%)')
plt.title('Variation within outcome due to Changes in Dose Variability',y=1.05)
plt.legend(loc='upper best')

