
# coding: utf-8

# ### Created all parts as functions to allow to be run multiple times with varied parameters

# ### Import Required Modules 

# In[1]:

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.optimize as opt
import pandas as pd

#%matplotlib qt
# qt if plot in seperate window
# inline if plot on page

## Set the number of d.p. displayed in numpy arrays
np.set_printoptions(precision=3)

## Progress bar widget
#from IPython.html.widgets import FloatProgress
from IPython.display import display

## allow timing of events
import time

## allow export of results as csv
import csv

## allow set of values to iterate over to be created
import itertools

## allow to get current workign directory for info.
import os

import datetime as dt

## Removed depreciation warning due to ints vs floats. Will not change results
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


# In[2]:

# ## For rounding to nearest nuber with specified precision use round_to(n,10) for nearest 10 for example.
# def round_to(n, precision):
#     correction = 0.5 if n >= 0 else -0.5
#     return int( n/precision+correction ) * precision

# def round_to_05(n):
#     return round_to(n, 0.05)


# In[3]:

### Alpha beta calculator - normal dist

def alphacalc_normal(alphabeta, sd):
    """Return alphabetanew and alpha from normal distribution as specified by sd.
    Default is beta = 0.03
    'alphabeta' is the alphabeta ratio
    If a negative value is returned it is resampled until positive"""
    
    beta = 0.03 # fixed beta in function
    
    ## get alpha beta to use from normal distribution
    if sd == 0:
        alphabetanew = alphabeta
    else:
        alphabetanew=np.random.normal(loc = alphabeta, scale = sd)
    
    ## make sure a positive value is returned
    while alphabetanew <= 0:
        alphabetanew=np.random.normal(loc = alphabeta, scale = sd)
    
    alpha = beta*alphabetanew
   
    return alpha, beta
## alpha/beta can be calced form the returned alpha and beta values


# In[4]:

### Alpha beta calculator - log normal dist

def alphacalc_lognormal(alphabeta, sd):
    """Return alphabetanew and alpha from normal distribution as specified by sd.
    Default is beta = 0.03
    'alphabeta' is the alphabeta ratio"""
    
    beta = 0.02 # fixed beta in function
    
    alphabeta_lognormal = np.log((alphabeta**2)/(np.sqrt((sd**2)+(alphabeta**2))))
    sd_lognormal = np.sqrt(np.log(((sd**2)/(alphabeta**2))+1))
    
    ## get alpha beta to use from normal distribution
    if sd == 0:
        alphabetanew = alphabeta
    else:
        alphabetanew=np.random.lognormal(mean = alphabeta_lognormal, sigma = sd_lognormal)
    
    alpha = beta*alphabetanew
   
    return alpha, beta
## alpha/beta can be calced form the returned alpha and beta values


# ###N0 varies depending on the a/b Std Dev!!!
# 
# Need to rethink how to determine N0 in a simple way.
# 
# - Ideally re-write entire thing to accept a set of parameters and then use these single parameters to feed into other functions.
# 
# 
# **Ideas**
# Could use some sort of fitting with all parameters fixed except N0.
# Then try and minimise the difference between TCP_calc and TCP_input?
# - Can set it to repeat the TCP calc for 1000 patients, then returnt he TCP at d_interest.
# - This would use just a single set of nominal parameters (as input by user/in model).
# - Only after this does the simulation of multiple parameter variations start.
# 
# I cannot think of another way of doing this, as I ahve been doing it manually (trial and error) so far...
# 

# In[5]:

## function to calculate the difference between input TCP and calculated TCP. only N0 is varied as want to estimate this.
def calc_dif_sq_orig(x,TCP, n, alphabeta_use, alphabeta_sd_use,d,d_shift,d_sd,d_trend,max_d,dose_of_interest,TCP_input):
    TCP_in = TCP
    
#    TCP_calc_all = completeTCPcalc(n=2000,
#                               alphabeta_use=3,
#                               alphabeta_sd_use=1,
#                               d=2,
#                               d_shift=0,
#                               d_sd=1,
#                               n0=x, # this needs to be able to vary to allow optimisation of the value
#                               max_d=100,
#                               dose_of_interest=74)
    TCP_calc_all = completeTCPcalc(n=n,
                               alphabeta_use=alphabeta_use,
                               alphabeta_sd_use=alphabeta_sd_use,
                               d=d,
                               d_shift=d_shift,
                               d_sd=d_sd,
                               d_trend=d_trend,
                               n0=x, # this needs to be able to vary to allow optimisation of the value
                               max_d=max_d,
                               dose_of_interest=dose_of_interest,
                               TCP_input=TCP_input)
    TCP_result = TCP_calc_all[10] ## Get only the result of interest (N0 is in posn. 10 int he returned tuple)
    
    TCP_Dif_sq = (TCP_in - TCP_result)**2 ## difference in squares to minimise
    return TCP_Dif_sq


# In[5]:

## function to calculate the difference between input TCP and calculated TCP. only N0 is varied as want to estimate this.
def calc_dif_sq(x,TCP, n, alphabeta_use, alphabeta_sd_use,d,d_shift,d_sd,d_trend,max_d,dose_of_interest,dose_input,TCP_input):
    
    ## tcp/dose_input must be provided as a list or this will not work (even if single valued).
    
    TCP_in = TCP_input
    #print(len(TCP_input))
    
    TCP_Dif_sq_sum = 0
    for i in range(len(TCP_input)):
        #print('TCP list val: ' + str(i))
        TCP_calc_all = completeTCPcalc(n=n,
                                   alphabeta_use=alphabeta_use,
                                   alphabeta_sd_use=alphabeta_sd_use,
                                   d=d,
                                   d_shift=d_shift,
                                   d_sd=d_sd,
                                   d_trend=d_trend,
                                   n0=x, # this needs to be able to vary to allow optimisation of the value
                                   max_d=max_d,
                                   dose_of_interest=dose_input[i],
                                   dose_input=dose_input[i],
                                   TCP_input=TCP_input[i])
        TCP_result = TCP_calc_all[10] ## Get only the result of interest (TCP at d_int is in posn. 10 in the returned tuple)
        TCP_Dif_sq = (TCP_in[i] - TCP_result)**2 ## difference in squares to minimise
        TCP_Dif_sq_sum = TCP_Dif_sq_sum + TCP_Dif_sq
        #print(TCP_Dif_sq_sum)
    return TCP_Dif_sq_sum

# In[6]:

## Determine N0 by minimising difference betwen calculate TCP and input TCP.

#def n0_determination_orig(TCP_input,
#                     repeats,
#                     n,
#                     alphabeta_use,
#                     alphabeta_sd_use,
#                     d,
#                     d_shift,
#                     d_sd,
#                     d_trend,
#                     max_d,
#                     dose_of_interest):
#    
#    ## TCP to fit N0 to. This will come from the literature.
#    
#    TCP_input = TCP_input
#    
#    ## Run optimisation multiple times to ensure accuracy (as based on random data)
#    repeats = repeats
#    n=n
#    alphabeta_use=alphabeta_use
#    alphabeta_sd_use=alphabeta_sd_use
#    d=d
#    d_shift=d_shift
#    d_sd=d_sd
#    d_trend=d_trend
#    max_d=max_d
#    dose_of_interest=dose_of_interest
#    
#    ## store the fit results in this list
#    fit_results=[]
#
#    #TCP, n, alphabeta_use, alphabeta_sd_use,d,d_shift,d_sd,max_d,dose_of_interest
#    
#    for i in range(0,repeats):
#        n0_result = opt.minimize_scalar(calc_dif_sq,method='brent',
#                                        args=(TCP_input,
#                                              n,
#                                              alphabeta_use,
#                                              alphabeta_sd_use,
#                                              d,
#                                              d_shift,
#                                              d_sd,
#                                              d_trend,
#                                              max_d,
#                                              dose_of_interest))
#        fit_results.append(n0_result.x) ## build up a list of the N0 fit results
#        no_fit_perc_completed = 100*i/repeats  ## Display the progress
#        print("\r" + "Fitting N0: " + str(int(no_fit_perc_completed)) + "% complete", end="")
#        
#    print('\r' + "Fitting N0: 100% complete")
#    #print("Fitting N0: 100% complete")
#
#    ## Remove outliers from N0 predictions (based on keeping the central 50% of results
#    trimmed = fit_results[:]
#    while len(trimmed)>0.5*len(fit_results):
#        trimmed.remove(max(trimmed))
#        trimmed.remove(min(trimmed))
#
#    ## Calculate mean of remaining results
#    n0_mean_fit=sum(trimmed)/len(trimmed)
#    
#    print("N0 mean fit: " + str(round(n0_mean_fit,2)))
#
#    ## Perform check calculation - not needed for actual use
#    TCP_Calc_min = completeTCPcalc(n=n,
#                               alphabeta_use=alphabeta_use,
#                               alphabeta_sd_use=alphabeta_sd_use,
#                               d=d,
#                               d_shift=d_shift,
#                               d_sd=d_sd,
#                               d_trend=d_trend,
#                               n0=n0_mean_fit, # this needs to be able to vary to allow optimisation
#                               max_d=max_d,
#                               dose_of_interest=dose_of_interest)
#    print('TCP target: ' + str(TCP_input))
#    print('TCP fit check: ' + str(round(TCP_Calc_min[10],2)))
#    
#    return n0_mean_fit, fit_results, trimmed
#
#
## ### a/b SD
## This parameter makes a significant difference to the N0 and the TCP values

# In[7]:

## N0 determination from entered TCP model parameters - DOES NOT WORK
#def n0_determination_old(TCP,d=2,D=74,ab=3,b=0.02):
#    """Determine the value of N0 from the input of the
#    nominal TCP model input parameters.
#    This enables the model to match that of clinical trials data if input."""
#    a = b * ab #ab ratio = a/b
#    n0 = -np.log(TCP)/np.exp(-(a*D)-(b*d*D))
#    return n0

#%%
def n0_determination(TCP_input,
                     n,
                     n0,
                     alphabeta_use,
                     alphabeta_sd_use,
                     d,
                     d_shift,
                     d_sd,
                     d_trend,
                     max_d,
                     dose_of_interest,
                     dose_input,
                     repeats = 5):

    ## TCP to fit N0 to. This will come from the literature.
    
## IF value of N0 is supplied, then do not need to fit it...
    #print('d_input: ' + str(dose_input))
    #print('dose_of_interest: ' +str(dose_of_interest))
   
    ## Run optimisation multiple times to ensure accuracy (as based on random data)
    repeats_min = 5
    if repeats < repeats_min:
        repeats = repeats_min
        print('Number of repeats for fitting N0 has been set to the minimum reccomended value of ' + str(repeats_min))
    else:
        repeats = repeats
    n=500 # set value of n for reliable fitting # use 3000 for final results
    n0=n0
    alphabeta_use=alphabeta_use
    alphabeta_sd_use=alphabeta_sd_use
    d=d
    d_shift=d_shift
    d_sd=d_sd
    d_trend=d_trend
    max_d=max_d
    dose_of_interest=dose_of_interest
    TCP_input = TCP_input
    dose_input = dose_input
    
    ## store the fit results in this list
    fit_results=[]

    #TCP, n, alphabeta_use, alphabeta_sd_use,d,d_shift,d_sd,max_d,dose_of_interest
    
### This is minimising the difference of squares returned by the function.
### This could probably be made to return the value if multiple TCP/Dose points are passed to it?
        
    print('Fitting N0 value')
    print('')
    
    for i in range(0,repeats):
        print('\r' + 'Fitting N0: Stage ' + str(i+1) + ' of ' + str(repeats), end="")
        n0_result = opt.minimize_scalar(calc_dif_sq,method='brent',
                                        args=(TCP_input,
                                        n,
                                        alphabeta_use,
                                        alphabeta_sd_use,
                                        d,
                                        d_shift,
                                        d_sd,
                                        d_trend,
                                        max_d,
                                        dose_of_interest,
                                        dose_input,
                                        TCP_input))
        #print(n0_result)
        fit_results.append(n0_result.x)
    fit_results.sort() # sort the repeated results to eliminate outliers
    #print(fit_results)
    #print(np.mean(fit_results))
    num_outliers = 1
    fit_results_trim = fit_results[num_outliers:-num_outliers]

    n0_mean_fit = sum(fit_results_trim)/len(fit_results_trim)
    #print(n0_mean_fit)
    print('')
    print('Fitting Completed')
    
    return n0_mean_fit#, TCP_Calc_min[10]

# In[8]:

# tester = np.array([])
# tester2 = np.array([])

# n2=10000

# mean_n = 10
# sd_n = 2

# ### mean and sd required for lognormal calc as described here: http://uk.mathworks.com/help/stats/lognstat.html?refresh=true
# mean_l = np.log((mean_n**2)/(np.sqrt((sd_n**2)+(mean_n**2))))
# sd_l = np.sqrt(np.log(((sd_n**2)/(mean_n**2))+1))

# for i in range(0,n2):
#     tester = np.append(tester,[np.random.normal(mean_n,sd_n)])
#     tester2 = np.append(tester2,[np.random.lognormal(mean_l,sd_l)])

# #tester = np.reshape(tester,(n2,12))
# #tester2 = np.reshape(tester2,(n2,2))
# #print(tester)

# print(sd_l)

# #plt.hist(tester, bins=30, color='red', alpha=0.5, label='norm')
# #plt.hist(tester2, bins=30, color='blue', alpha=0.5, label='lognorm')
# #plt.legend()
# ##print(tester2)
# #plt.show()


# In[9]:

###calculate dose for a given fraction based on normal distribution around dose shift
### Log normal distribution not necessary as SD small wrt mean so v unlikely to get negative values returned.
### Could also set a limit on the range (say to +/- 5%)? But wont make much difference with a small SD anyway.

def fracdose(dose, shift, sd):
    """Return dose_actual from normal distribution around dose (Gy) as specified by sd (%) and shift (%).
    Default is dose = 2Gy, shift = 0%, and sd of 0%
    If a negative value is returned it is resampled until positive
    The standard deviation is of the nominal dose"""
    
    ## get actual dose to use from normal distribution based on shift
    
    dose_shift = dose + (dose*shift/100)
    
    ## if sd is zero, then no change to dose
    if sd == 0:
        dose_actual = dose_shift
        return dose_actual
    
    dose_actual=np.random.normal(loc = dose_shift, scale = (dose*sd/100))
    
    ## make sure a positive value is returned
    while dose_actual <= 0:
        dose_actual=np.random.normal(loc = dose_shift, scale = (dose*sd/100))
    
    return dose_actual


# In[10]:

## Survival Fraction Calculation
def SFcalc(alpha, beta, dose):
    """Return the SF with input values.
    Note this is for a single dose delivery.
    The product of multiple fractions shoudld be taken
    to give overall SF"""
    
    SF = np.exp(-(alpha*dose) - (beta*(dose**2)))
    
    return SF


# In[11]:

## TCP Calculation absed on cumulative SF
def TCPcalc(sf, n0):
    """Return the TCP with input values.
    Based on cumulative SF and N0"""
    
    TCP = np.exp(-n0*sf)
    
    return TCP


# In[12]:

## Calc Number of fractions to get to max dose (note: round up as always want an integer)

def no_frac_nom_doses_array(max_d, d):
    n_frac = np.ceil(max_d/d)

    fractions = np.arange(1,n_frac+1)
    #print(fractions)
    nom_doses = np.arange(d,(d*n_frac)+d, step = d)
    #print(nom_doses)
    return fractions, nom_doses, n_frac


# In[13]:

## This gives a column with the patient number and makes it easier to check values as going

def create_patients(n):
    
    if n<1:
        n=1
    patients = np.arange(0,n)+1
    patients.shape=(n,1)
    #print(patients)
    return patients


# In[14]:

## empty array to store alpha values in (log normal distribution to be used)

def create_alpha_beta_array(n, alphabeta_use, alphabeta_sd_use):
    alpha_and_beta =[]

    for p in range(0,n):
        #alpha_and_beta = np.append(alpha_and_beta,[alphacalc_normal(alphabeta = alphabeta_use, sd=alphabeta_sd_use)])
        alpha_and_beta.append([alphacalc_lognormal(alphabeta = alphabeta_use, sd=alphabeta_sd_use)])

    ## reshape to get a row per patient
    alpha_and_beta_np = np.array(alpha_and_beta)
    alpha_and_beta_np = np.reshape(alpha_and_beta_np,(n,2))
    #print(alpha_and_beta)
    return alpha_and_beta_np


# In[15]:

## Calculate Doses for all patients and all fractions and put in an array
## Added in a linear output trend. Specify in percent per day/treatment. Default = 0

def doses_array(n, n_frac, d, d_shift, d_sd, d_trend=0):
    doses = []
    for i in range(0,int(n*n_frac)):
        doses.append(fracdose(dose = d, shift=d_shift, sd=d_sd))
    doses_np = np.array(doses)
    doses_np = np.reshape(doses_np,(n,n_frac))
    
    ## Make an array of the trend to multiply the doses array with
    trend_array = []
    for i in range(int(n_frac)):
        trend_array.append(1+(i*d_trend/100))
        
    ## multiply array of doses by trend value for each fraction    
    doses_np_trend = doses_np*trend_array
        
    #print(doses)
    return doses_np_trend
    

# In[16]:

# ## Example calc and plot of doses with trend.

# x = doses_array(10,100,2,-2,0.5,3/365)
# print(x[1])
# for i in range(9):
#     plt.plot(x[i], alpha=0.5)
    
# plt.plot(x[1], linewidth=3, c='black')
# plt.show()


# In[17]:

## Calculate Doses for all patients and all fractions and put in an array

#def doses_array_orig(n, n_frac, d, d_shift, d_sd):
#    doses = []
#    for i in range(0,int(n*n_frac)):
#        doses.append(fracdose(dose = d, shift=d_shift, sd=d_sd))
#    doses_np = np.array(doses)
#    doses_np = np.reshape(doses_np,(n,n_frac))
#    #print(doses)
#    return doses_np


# In[18]:

## Combine all results into single array which may be easier to work with for analysis

def combine_results(patients, alpha_and_beta, doses):
    results_wanted = (patients,alpha_and_beta, doses)
    all_results = np.concatenate(results_wanted, axis=1)
    #print(all_results)
    return all_results


# In[19]:

## Loop through the doses of the first patient (first row [0] of array)
def calc_all_SFs(patients, n, n_frac, alpha_and_beta, doses):
    SFs = []

    for i in range(0,len(patients)): # loop through each patient (row)
        for j in range(0,int(n_frac)): # loop through each fraction for each patient (col)
            SFs.append(SFcalc(alpha_and_beta[i][0],alpha_and_beta[i][1],doses[i,j]))

    SFs_np = np.array(SFs)
    SFs_np = np.reshape(SFs_np,(n,n_frac))

    ## GEt cumulative SF for each patient
    SF_cum = np.cumprod(SFs_np, axis=1)
    return SFs_np, SF_cum


# In[20]:

## append results to text file.

def saveasCSV(filename, array):

    fl = open(filename, 'a', newline='\n')

    writer = csv.writer(fl)
    writer.writerow(array)

    fl.close()


#%%
def d_list_sort(d_list,n_frac,n):
    d_list_pad = d_list + [0] * (n_frac - len(d_list))# pad with zero doses 
    #print(d_list_pad)
    doses = [d_list_pad for i in range(n)] # set the provided list as the doses
    doses = np.array(doses) # convert to numpy array for use to with other data types
    return doses


# In[21]:

## Calc Number of fractions and nominal dose per fraction to get to max dose
    
def completeTCPcalc(n,
                   alphabeta_use,
                   alphabeta_sd_use,
                   d,
                   d_shift,
                   d_sd,
                   d_trend,
                   max_d,
                   dose_of_interest,
                   dose_input = None,
                   TCP_input=None,
                   d_list=None,
                   n0=None):
    #print(dose_input)
    fractions, nom_doses, n_frac = no_frac_nom_doses_array(max_d,d)

    if d_list is not None:
        doses = d_list_sort(d_list,n_frac,n)
        #d_list_pad = d_list + [0] * (n_frac - len(d_list))# pad with zero doses 
        #print(d_list_pad)
        #doses = [d_list_pad for i in range(n)] # set the provided list as the doses
        #doses = np.array(doses) # convert to numpy array for use to with other data types
    else:
        doses = doses_array(n, n_frac, d, d_shift, d_sd, d_trend)
    #print(doses)
    #**
    ## array of doses after each fraction for each patient
        
    ## create array containing number of patients in population
    patients = create_patients(n)
    #print(patients)
    
    ## Creat array of alpha and veta values for each patient
    alpha_and_beta = create_alpha_beta_array(n, alphabeta_use, alphabeta_sd_use)
    
    ## put all results in an array with a patient on each row
    all_results = combine_results(patients, alpha_and_beta, doses)

    ## Calc cumulative SF for all patients (also return individual fraction SFs)
    SFs, SF_cum = calc_all_SFs(patients, n, n_frac, alpha_and_beta, doses)
    
    ## determine N0 value to use if none provided
    if (n0 is None) and (TCP_input is not None) and (dose_input is not None):
        print('N0 not provided')
        n0_use = n0_determination(TCP_input,
                                  n,
                                  n0,
                                  alphabeta_use,
                                  alphabeta_sd_use,
                                  d,
                                  d_shift,
                                  d_sd,
                                  d_trend,
                                  max_d,
                                  dose_of_interest,
                                  dose_input)
    else:
        n0_use = n0

    ## Calculate TCP for all individual patients and fractions
    TCPs = TCPcalc(sf = SF_cum, n0=n0_use)

    ## Calculate population TCP by averaging down the columns
    TCP_pop = np.mean(TCPs, axis = 0)

    frac_of_interest = dose_of_interest/d #reinstate this
    #frac_of_interest = 74/d

    TCP_at_dose_of_interest = TCP_pop[frac_of_interest]

    #t_end = time.time()

    #t_total = t_end-t_start
    #print(str(round(t_total,2)) + ' secs for ' + str(n) + ' patients')

    TCPs_of_interest = TCPs[:,frac_of_interest-1]

    TCP_cure = (TCPs_of_interest).sum()
    TCP_cure_percent = 100*TCP_cure/n
    #print(TCP_cure_percent)
    #print(TCP_pop)
    
    return n,alphabeta_use,alphabeta_sd_use,d,d_shift,d_sd,n0_use,max_d,dose_of_interest,frac_of_interest,TCP_cure_percent, TCPs, TCP_pop, nom_doses, d_trend


# In[22]:

## Plot of individual and population TCPs as a function for ease

def TCP_plot(no_ind_plots, label):
    #no_ind_plots = 50

    ## individual plots cannot be more than total patients
    if(no_ind_plots>n):
        no_ind_plots=n

    ## want to select the individual plots randomly from those calcualted...
    ind_plots = np.random.choice(len(TCPs),no_ind_plots, replace=False)

    ## individuals (specified number of plots chosen)
    for i in ind_plots:
        plt.plot(nom_doses,TCPs[i], color = 'grey', alpha = 0.5)
    ## population
    plt.plot(nom_doses,TCP_pop, color='black', linewidth='2', alpha=0.5)
    plt.plot(nom_doses,TCP_pop, marker = 'o', ls='none', label=label)

    ## plot formatting
    plt.xlim(0,max(nom_doses))
    plt.ylim(0,1.0)
    plt.xlabel('Dose (Gy)')
    plt.ylabel('TCP')
    plt.title('TCPs')
    #plt.legend(loc = 'best', fontsize = 'medium', framealpha = 1)
    plt.axvline(d_interest, color = 'black', ls='--',)
    plt.axhline(TCP_pop[frac_interest-1], color='black', ls='--')

    ## add labels with TCP at dose of interest
    text_string = ('Pop. TCP = ' + str(round(TCP_cure_at_d_interest,2)) + ' % at ' + str(d_interest) + 'Gy')
    plt.text(5,0.4,text_string, backgroundcolor='white')
    plt.legend(loc = 'lower left',numpoints=1)

    plt.show()


# In[23]:

## Return sum of squares [calc = calculation of TCP from TCP calc, ref = TCP_input]
def sq_dif(calc, ref):
    sq_dif = (calc - ref)**2
    return sq_dif


# In[24]:

# ## Determine N0 by minimising difference betwen calculated TCP and input TCP.

# TCP_input = 80 # TCP in percent

# ### Need to add something to ensure the minimum is always found, but without always having a huge list of possible values.
# ## Could also start with a coarse list and refine after one sweep?
# j=1
# start = 0
# increment=10

# any_pos = False
# any_neg = False
# stop = False

# while (stop == False and j<100):
#     TCP_list = []
       
#     a_list = list(range(start,(j*increment)))
#     print(a_list)
    
#     for i in a_list:
#         TCP_Calc = completeTCPcalc(n=100,
#                                alphabeta_use=3,
#                                alphabeta_sd_use=1,
#                                d=2,
#                                d_shift=0,
#                                d_sd=1,
#                                d_trend=0,
#                                n0=i, # this needs to be able to vary to allow optimisation
#                                max_d=100,
#                                dose_of_interest=74)
#         TCP_list.append(TCP_Calc[10])

#     dif_list = np.array(TCP_list)-TCP_input
#     dif_list_abs = abs(np.array(TCP_list)-TCP_input)
#     print("test")
#     print(dif_list)

#     any_pos = np.any(dif_list >0)
#     any_neg = np.any(dif_list <0)
    
#     if (any_pos == True and any_neg == True):
#         stop = True
#         #break
    
#     j = j+1 # increment list to lookup
#     start = start + increment
    
# print("N0 fitted")
# ## can use the following to check if the minimum has been found by cehckign for + and - differences.


# #dif_list[:] = [x-TCP_input for x in TCP_list]

# #a[:] = [x - 13 for x in a]
    
# ### Now want to 

# #tester = sp.optimize.minimize(sq_dif(TCP_Calc(x),TCP_input))
# print(min(dif_list_abs)) # minimum difference found.
# print(np.argmin(dif_list_abs))# index of minimum found
# print("N0: " + str(a_list[np.argmin(dif_list_abs)]))

# ## So n0 = a_list[np.argmin(dif_list)]

# TCP_Calc_min = completeTCPcalc(n=1000,
#                            alphabeta_use=3,
#                            alphabeta_sd_use=1,
#                            d=2,
#                            d_shift=0,
#                            d_sd=1,
#                            d_trend=0,
#                            n0=a_list[np.argmin(dif_list_abs)], # this needs to be able to vary to allow optimisation
#                            max_d=100,
#                            dose_of_interest=74)
# print(TCP_Calc_min[10])


# In[25]:

# n0_mean_fit, fit_results, trimmed = n0_determination(TCP_input=80,
#                                                      repeats=20,
#                                                      n=1000,
#                                                      alphabeta_use=3,
#                                                      alphabeta_sd_use=1,
#                                                      d=2,
#                                                      d_shift=0,
#                                                      d_sd=1,
#                                                      d_trend=0,
#                                                      max_d=100,
#                                                      dose_of_interest=74)


# In[26]:

# ## Plot of N0 fit results
# plt.hist(fit_results)
# plt.hist(trimmed)
# plt.show()


# ### Calculate and Plot results of single set of parameters

# In[27]:

# ### Individual calculation of TCP at a set value of dose

# n_rpts = 1

# print("Enter Dose Shift (%)")
# applied_shift = float(input())

# text_string = ("Shift: " + str(applied_shift) + "%")

# for i in range(0,n_rpts):
#     t = completeTCPcalc(n = 1000,
#                     alphabeta_use = 3,
#                     alphabeta_sd_use = 2,
#                     d = 2,
#                     d_shift = applied_shift, # +1.6,-2.1% from OP data
#                     d_sd = 0,
#                     d_trend=0,
#                     n0 = n0_mean_fit,
#                     max_d = 100,
#                     dose_of_interest = 74)

#     n=t[0]
#     TCP_pop = t[-2]
#     TCPs = t[-3]
#     nom_doses = t[-1]
#     d_interest = t[8]
#     frac_interest = t[9]
#     TCP_cure_at_d_interest = t[10]



#     #print(TCP_pop)
#     TCP_plot(100,text_string)


# ### User iterator to build list to vary multiple parameters in one go

# In[28]:

## Round to nearest n (e.g. nearest 50 n=50)
def round_n(x,n):
    return int(round(x / n)) * n


# In[29]:

######***** Change this to calcualte a given percentage variation in dose? or allow both.....?

## This allows simple building of a tuple containing all possible combinations of values

def dose_iter(dose_var=0.5,
             dose_max=3,
             ab_var=1,
             ab_max=4,
             ab_min=2,
             ab_sd_var=0.5,
             ab_sd_max=1,
             ab_sd_min=0.5,
             di_var=0,
             di_max=74,
             di_min=74,
             n0_nominal=160,
             n0_var=5,
             n0_range=10):
    
    ## dose shifts to test
    dose_var = dose_var # dose step size
    dose_max = dose_max
    dose_min = -dose_max
    dose_number = (dose_max-dose_min)/dose_var+1 # number of points

    dose_vals = np.linspace(dose_min,dose_max,dose_number) # dose values to test

    ## alphabeta values to test
    ab_var = ab_var # step size
    ab_max = ab_max
    ab_min = ab_min
    ab_number = (ab_max-ab_min)/ab_var+1 # number of points
    ab_vals = np.linspace(ab_min,ab_max,ab_number) # different a/b values to test

    ## alphabeta_SD values to test
    ab_sd_var = ab_sd_var # step size
    ab_sd_max = ab_sd_max
    ab_sd_min = ab_sd_min
    ab_sd_number = (ab_sd_max-ab_sd_min)/ab_sd_var+1 # number of points
    ab_sd_vals = np.linspace(ab_sd_min,ab_sd_max,ab_sd_number) # different a/b values to test


    ## dose of interest to test
    di_var = di_var # dose step size
    di_max = di_max
    di_min = di_min
    di_number = (di_max-di_min)/di_var+1 # number of points
    di_vals = np.linspace(di_min,di_max,di_number)

    ## N0 values - determined from TCP model paramters input (from clinical data - coudl get form text file etc?)
    #n0_nominal = n0_determination(TCP=0.88,d=2,D=74,ab=3,b=0.02)
    n0_nominal = n0_nominal
    n0_var = n0_var
    n0_range = n0_range #percentage variation
    n0_min = round_n(n0_nominal * (1-n0_range/100),n0_var)
    n0_max = round_n(n0_nominal * (1+n0_range/100),n0_var)
    n0_number = (n0_max-n0_min)/n0_var+1
    n0_vals = np.linspace(n0_min,n0_max,n0_number)

    total_tests = len(dose_vals)*len(ab_vals)*len(di_vals)*len(n0_vals)*len(ab_sd_vals)
    test_val_iterator = itertools.product(dose_vals,ab_vals,di_vals,n0_vals,ab_sd_vals)

    test_vals = list(test_val_iterator) #list of all combinations of parameters
    num_its = len(test_vals) #number of iterations
    
    return test_vals,num_its
    
#print(num_its)
#print(n0_mean_fit)
#print(n0_vals)


# In[30]:

# abc1 = dose_iter(dose_var=0.5,
#              dose_max=3,
#              ab_var=0.5,
#              ab_max=4,
#              ab_min=2,
#              di_var=1,
#              di_max=80,
#              di_min=64,
#              n0_var=10,
#              n0_range=5)

# print(abc1[1])


# In[38]:

## vary multiple values through use of constructer iterator list above.

#### This needs the N0 to be calculated before calculation starts.
### The values for the iterator above shouls also be calculated within this?
#t1 = datetime.now()

def TCP_full(k=10,
             TCP_input=80,
             repeats=20,
             n=1000,
             n0=169,
             alphabeta_use=3,
             alphabeta_sd_use=1,
             d=2,
             d_shift=0,
             d_sd=1,
             d_trend=0,
             max_d=100,
             dose_of_interest=74,
             save_name="TCPresultsProst-all_results-std.csv"):
    
    # gat variable values (does this really need doing???)
    TCP_input=TCP_input
    repeats=repeats
    n=n
    n0=n0
    alphabeta_use=alphabeta_use
    alphabeta_sd_use=alphabeta_sd_use
    d=d
    d_shift=d_shift
    d_sd=d_sd
    d_trend=d_trend
    max_d=max_d
    dose_of_interest=dose_of_interest
    save_name=save_name+".csv"
    
    #Print some values
    print("TCP Point at dose of interest: " + str(TCP_input))
    print("Dose of interest: " + str(dose_of_interest))
    
    print("Starting TCP Simulation")
    
    #N0 determination should use nominal values of dose to assume average?
    ###orig has n0_mean_fit, fit_results, trimmed = n0_determ....
    if n0 is None:
        n0_mean_fit = n0_determination(TCP_input=TCP_input,
                                       repeats=repeats,
                                       n=n,
                                       n0=n0,
                                       alphabeta_use=alphabeta_use,
                                       alphabeta_sd_use=alphabeta_sd_use,
                                       d=d,
                                       d_shift=0,
                                       d_sd=0,
                                       d_trend=0,
                                       max_d=max_d,
                                       dose_of_interest=dose_of_interest)
    else:
        n0_mean_fit = n0
    
    set_precision = 5 #round to nearest 5
    n0_mean_fit_set_precision = round_n(n0_mean_fit,set_precision)
    
    ab_range = 1
    
    iter_list = dose_iter(dose_var=0.5,
                          dose_max=4,
                          ab_var=1,
                          ab_max=alphabeta_use + ab_range,
                          ab_min=alphabeta_use - ab_range,
                          di_var=2,
                          di_max=80,
                          di_min=64,
                          n0_nominal=n0_mean_fit_set_precision,
                          n0_var=5,
                          n0_range=10)
    
    print("N0 mean fit (nearest " + str(set_precision) + "): " + str(n0_mean_fit_set_precision))
    
    test_vals = iter_list[0]
    num_its = len(test_vals)
    print("Total simulations to run: " + str(num_its))
    
    k = k # number of repeats

    #f = FloatProgress(min=0, max=num_its)
    #display(f)

    #f1 = FloatProgress(min=0, max=k-1)
    #display(f1)

    #barpos = 1

    all_test_results_array = np.array([]) # array to store all results in before saving as excel csv file
    all_test_results_array = []
    
    #print(test_vals)
    start_time = dt.datetime.now()
    progress_perc = 0
    for j in test_vals:

        current_time = dt.datetime.now()

        progress_perc = progress_perc + (100/len(test_vals))
        no_completed = progress_perc/100 * len(test_vals)
        
        dif_secs = (current_time-start_time).seconds
        
        remaining_secs = (num_its-no_completed)*(dif_secs/no_completed)
        remaining_mins = remaining_secs/60
        #print(remaining_secs)
        
        #print('\r' + "TCP Calculation: " + str(int(progress_perc)) + "% completed", end='')
        print('\r' + "Running TCP Simulation: " + str(round(no_completed,0)) + " of " + str(len(test_vals)) + " (" + str(round(progress_perc,1)) + "% completed)" + " (" + str(round(remaining_mins,1)) + " mins remaining)", end='')

        results_array = []

        n = n
        
        for i in range(0,k):
            t = completeTCPcalc(n,
                            alphabeta_use = j[1],
                            alphabeta_sd_use = j[4],
                            d = d,
                            d_shift = j[0],
                            d_sd = d_sd,
                            d_trend=d_trend,
                            n0 = j[3],
                            max_d = max_d,
                            dose_of_interest = j[2])
            
            #f1.value = i

            n=t[0]
            alphabeta_use = t[1]
            alphabeta_sd_use = t[2]
            d = t[3]
            d_shift = t[4]
            d_sd = t[5]
            d_trend = t[-1]
            n0 = t[6]
            TCP_pop = t[-3]
            TCPs = t[-4]
            nom_doses = t[-2]
            d_interest = t[8]
            frac_interest = t[9]
            TCP_cure_at_d_interest = t[10]
            max_d = max_d

            results_array.append(TCP_cure_at_d_interest)

            param_array = []
            param_array.append(n)
            param_array.append(k)
            param_array.append(alphabeta_use)
            param_array.append(alphabeta_sd_use)
            param_array.append(d)
            param_array.append(d_shift)
            param_array.append(d_sd)
            param_array.append(d_trend)
            param_array.append(n0)
            param_array.append(max_d)
            param_array.append(d_interest)
            
            header_array = []
            header_array.append("n")
            header_array.append("k")
            header_array.append("ab")
            header_array.append("ab_sd")
            header_array.append("d")
            header_array.append("d_shift")
            header_array.append("d_sd")
            header_array.append("d_trend")
            header_array.append("n0")
            header_array.append("max_d")
            header_array.append("d_interest")
            for i in range(k):
                header_array.append("k_"+str(i+1))
            

        ## Array containing results form a single set of parameters (this will update each iteration)
        param_array.extend(results_array)
        param_results_array = param_array

        ## Create an array containing all the results for all sets of parameters
        all_test_results_array.extend(param_results_array)

        # for updating progress bar
        #barpos = barpos+1
        #f.value = barpos

    #print(param_results_array[-1])

    ## Covnert to Numpy array and reshape so a single set of parameters is on one row
    all_test_results_array_np = np.array(all_test_results_array)
    all_test_results_array_np = np.reshape(all_test_results_array_np,(num_its,len(all_test_results_array_np)/num_its))

    ## Save results to file
    #np.savetxt(save_name, all_test_results_array_np, delimiter=",",fmt='%10.3f')
    
    with open(save_name, "w") as f:
        writer = csv.writer(f)
        writer.writerows([header_array])
        writer.writerows(all_test_results_array_np)
  
    print('\r' + "TCP Simulation: 100% Completed")
    
    print("Results saved in file: " + os.getcwd()+ save_name)
    #print(header_array)
    
    return save_name, n0_mean_fit

#t2 = datetime.now()
#t3 = (t2-t1).total_seconds()

#print(t3)


# In[39]:

## This runs the entire analysis in one go.
#all_test = TCP_full(k=5,
#                    TCP_input=80,
#                    repeats=5,
#                    n=100,
#                    alphabeta_use=3,
#                    alphabeta_sd_use=0.5,
#                    d=2,
#                    d_shift=0,
#                    d_sd=0.8,
#                    d_trend=0,
#                    max_d=100,
#                    dose_of_interest=74,
#                    save_name="Results1-Trend=0")


# ### Import created results to examine

# In[42]:

#TCP_data = pd.read_csv(os.getcwd()+"\\"+all_test[0]) # read in file saved after TCP simulation.


# In[43]:

#print(all_test[0])
#TCP_data.describe()


# In[ ]:



