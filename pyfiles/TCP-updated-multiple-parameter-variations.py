
# coding: utf-8

# ### Created all parts as functions to allow to be run multiple times with varied parameters

# ### Import Required Modules 

# In[1]:

import numpy as np
import matplotlib.pyplot as plt

get_ipython().magic('matplotlib qt')
# qt if plot in seperate window
# inline if plot on page

## Set the number of d.p. displayed in numpy arrays
np.set_printoptions(precision=3)

## Progress bar widget
from IPython.html.widgets import FloatProgress
from IPython.display import display

## allow timing of events
import time

## allow export of results as csv
import csv


# In[2]:

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


# In[3]:

### Alpha beta calculator - log normal dist

def alphacalc_lognormal(alphabeta, sd):
    """Return alphabetanew and alpha from normal distribution as specified by sd.
    Default is beta = 0.03
    'alphabeta' is the alphabeta ratio"""
    
    beta = 0.03 # fixed beta in function
    
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


# In[4]:

tester = np.array([])
tester2 = np.array([])

n2=10000

mean_n = 10
sd_n = 2

### mean and sd required for lognormal calc as described here: http://uk.mathworks.com/help/stats/lognstat.html?refresh=true
mean_l = np.log((mean_n**2)/(np.sqrt((sd_n**2)+(mean_n**2))))
sd_l = np.sqrt(np.log(((sd_n**2)/(mean_n**2))+1))

for i in range(0,n2):
    tester = np.append(tester,[np.random.normal(mean_n,sd_n)])
    tester2 = np.append(tester2,[np.random.lognormal(mean_l,sd_l)])
    
#tester = np.reshape(tester,(n2,12))
#tester2 = np.reshape(tester2,(n2,2))
#print(tester)

print(sd_l)

plt.hist(tester, bins=30, color='red', alpha=0.5, label='norm')
plt.hist(tester2, bins=30, color='blue', alpha=0.5, label='lognorm')
plt.legend()
#print(tester2)
plt.show()


# In[5]:

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


# In[6]:

## Survival Fraction Calculation
def SFcalc(alpha, beta, dose):
    """Return the SF with input values.
    Note this is for a single dose delivery.
    The product of multiple fractions shoudld be taken
    to give overall SF"""
    
    SF = np.exp(-(alpha*dose) - (beta*(dose**2)))
    
    return SF


# In[7]:

## TCP Calculation absed on cumulative SF
def TCPcalc(sf, n0):
    """Return the TCP with input values.
    Based on cumulative SF and N0"""
    
    TCP = np.exp(-n0*sf)
    
    return TCP


# In[8]:

## Calc Number of fractions to get to max dose (note: round up as always want an integer)

def no_frac_nom_doses_array(max_d, d):
    n_frac = np.ceil(max_d/d)

    fractions = np.arange(1,n_frac+1)
    #print(fractions)
    nom_doses = np.arange(d,(d*n_frac)+d, step = d)
    #print(nom_doses)
    return fractions, nom_doses, n_frac


# In[9]:

## This gives a column with the patient number and makes it easier to check values as going

def create_patients(n):
    
    if n<1:
        n=1
    patients = np.arange(0,n)+1
    patients.shape=(n,1)
    #print(patients)
    return patients


# In[10]:

## empty array to store alpha values in (log normal distribution to be used)

def create_alpha_beta_array(n, alphabeta_use, alphabeta_sd_use):
    alpha_and_beta = np.array([])

    for p in range(0,n):
        #alpha_and_beta = np.append(alpha_and_beta,[alphacalc_normal(alphabeta = alphabeta_use, sd=alphabeta_sd_use)])
        alpha_and_beta = np.append(alpha_and_beta,[alphacalc_lognormal(alphabeta = alphabeta_use, sd=alphabeta_sd_use)])

    ## reshape to get a row per patient
    alpha_and_beta = np.reshape(alpha_and_beta,(n,2))
    #print(alpha_and_beta)
    return alpha_and_beta


# In[11]:

## Calculate Doses for all patients and all fractions and put in an array

def doses_array(n, n_frac, d, d_shift, d_sd):
    doses = np.array([])
    for i in range(0,int(n*n_frac)):
        doses = np.append(doses,fracdose(dose = d, shift=d_shift, sd=d_sd))
    doses = np.reshape(doses,(n,n_frac))
    #print(doses)
    return doses


# In[12]:

## Combine all results into single array which may be easier to work with for analysis

def combine_results(patients, alpha_and_beta, doses):
    results_wanted = (patients,alpha_and_beta, doses)
    all_results = np.concatenate(results_wanted, axis=1)
    #print(all_results)
    return all_results


# In[13]:

## Loop through the doses of the first patient (first row [0] of array)

def calc_all_SFs(patients, n, n_frac, alpha_and_beta, doses):
    SFs = np.array([])

    for i in range(0,len(patients)): # loop through each patient (row)
        for j in range(0,int(n_frac)): # loop through each fraction for each patient (col)
            SFs = np.append(SFs,SFcalc(alpha_and_beta[i][0],alpha_and_beta[i][1],doses[i,j]))

    SFs = np.reshape(SFs,(n,n_frac))

    ## GEt cumulative SF for each patient
    SF_cum = np.cumprod(SFs, axis=1)
    return SFs, SF_cum


# In[14]:

## append results to text file.

def saveasCSV(filename, array):

    fl = open(filename, 'a', newline='\n')

    writer = csv.writer(fl)
    writer.writerow(array)

    fl.close()


# In[15]:

## Calc Number of fractions and nominal dose per fraction to get to max dose
    
def completeTCPcalc(n,
                    alphabeta_use,
                    alphabeta_sd_use,
                    d,
                    d_shift,
                    d_sd,
                    n0,
                    max_d,
                    dose_of_interest):

    fractions, nom_doses, n_frac = no_frac_nom_doses_array(max_d, d)

    ## create array containing number of patients in population
    patients = create_patients(n)

    ## Creat array of alpha and veta values for each patient
    alpha_and_beta = create_alpha_beta_array(n, alphabeta_use, alphabeta_sd_use)

    ## array of doses after each fraction for each patient
    doses = doses_array(n, n_frac, d, d_shift, d_sd)

    ## put all results in an array with a patient on each row
    all_results = combine_results(patients, alpha_and_beta, doses)

    ## Calc cumulative SF for all patients (also return individual fraction SFs)
    SFs, SF_cum = calc_all_SFs(patients, n, n_frac, alpha_and_beta, doses)

    ## Calculate TCP for all individual patients and fractions
    TCPs = TCPcalc(sf = SF_cum, n0=n0)

    ## Calculate population TCP by averaging down the columns
    TCP_pop = np.mean(TCPs, axis = 0)

    frac_of_interest = dose_of_interest/d

    TCP_at_dose_of_interest = TCP_pop[frac_of_interest]

    #t_end = time.time()

    #t_total = t_end-t_start
    #print(str(round(t_total,2)) + ' secs for ' + str(n) + ' patients')

    TCPs_of_interest = TCPs[:,frac_of_interest-1]

    TCP_cure = (TCPs_of_interest).sum()
    TCP_cure_percent = 100*TCP_cure/n
    
    return n,alphabeta_use,alphabeta_sd_use,d,d_shift,d_sd,n0,max_d,dose_of_interest,frac_of_interest,TCP_cure_percent, TCPs, TCP_pop, nom_doses


# In[55]:

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


# ### Calcualte and Plot results of single set of parameters

# In[88]:

## Individual calculation of TCP at a set value of dose

n_rpts = 1

for i in range(0,n_rpts):
    t = completeTCPcalc(n = 500,
                    alphabeta_use = 6,
                    alphabeta_sd_use = 2,
                    d = 3,
                    d_shift = 0,
                    d_sd = 0,
                    n0 = 1E5,
                    max_d = 100,
                    dose_of_interest = 60)

    n=t[0]
    TCP_pop = t[-2]
    TCPs = t[-3]
    nom_doses = t[-1]
    d_interest = t[8]
    frac_interest = t[9]
    TCP_cure_at_d_interest = t[10]
   
    #print(TCP_pop)
    TCP_plot(50,"Standard")


# ### Repeat calculation a set number of times

# In[16]:

## Multiple calcualtions of TCP with single set values

k = 20 # number of repeats

results_array = np.array([])

f = FloatProgress(min=0, max=k-1)
display(f)

for i in range(0,k):

    t = completeTCPcalc(n = 200,
                    alphabeta_use = 10,
                    alphabeta_sd_use = 2,
                    d = 2,
                    d_shift = 5,
                    d_sd = 0,
                    n0 = 1E9,
                    max_d = 100,
                    dose_of_interest = 74)
    
    n=t[0]
    alphabeta_use = t[1]
    alphabeta_sd_use = t[2]
    d = t[3]
    d_shift = t[4]
    d_sd = t[5]
    n0 = t[6]
    TCP_pop = t[-2]
    TCPs = t[-3]
    nom_doses = t[-1]
    d_interest = t[8]
    frac_interest = t[9]
    TCP_cure_at_d_interest = t[10]
    max_d = 100
    
    results_array = np.append(results_array,TCP_cure_at_d_interest)
    f.value = i # for updating progress bar
    
## include model parameters into array and insert before results
param_array = np.array([])
param_array = np.append(param_array,n)
param_array = np.append(param_array,k)
param_array = np.append(param_array,alphabeta_use)
param_array = np.append(param_array,alphabeta_sd_use)
param_array = np.append(param_array,d)
param_array = np.append(param_array,d_shift)
param_array = np.append(param_array,d_sd)
param_array = np.append(param_array,n0)
param_array = np.append(param_array,max_d)
param_array = np.append(param_array,d_interest)
#param_array = np.append(param_array,frac_of_interest)

param_results_array = np.concatenate((param_array,results_array))

#print(results_array)
#print(param_array)
print(param_results_array[-1])

#saveasCSV(filename = 'TCPresults.csv', array = param_results_array)


# ### User iterator to build list to vary multiple parameters in one go

# In[19]:

## will need to iterate over all the values and add them to a list.

## This allows simple building of a tuple containing all possible combinations of values
import itertools

## dose shifts to test
dose_var = 1 # dose step size
dose_max = 5
dose_min = -dose_max
dose_number = (dose_max-dose_min)/dose_var+1 # number of points

dose_vals = np.linspace(dose_min,dose_max,dose_number)
#print(dose_vals)

## alphabeta values to test
ab_var = 2 # dose step size
ab_max = 12
ab_min = 8
ab_number = (ab_max-ab_min)/ab_var+1 # number of points

ab_vals = np.linspace(ab_min,ab_max,ab_number)
#print(ab_vals)

## dose of interest to test
di_var = 4 # dose step size
di_max = 80
di_min = 70
di_number = (di_max-di_min)/di_var+1 # number of points

di_vals = np.linspace(di_min,di_max,di_number)
#print(di_vals)

total_tests = len(dose_vals)*len(ab_vals)*(len(di_vals))
#print(total_tests)

test_val_iterator = itertools.product(dose_vals,ab_vals,di_vals)

test_vals = list(test_val_iterator)

num_its = len(test_vals)

#for i in test_vals:
#    print(i)
    
print(num_its)


# In[20]:

## vary multiple values through use of constructer iterator list above.

k = 20 # number of repeats

#test_vals = ((-3,10),(-2,8),(-1,8),(0,10),(1,8),(2,8),(3,10)) # values to vary and test

f = FloatProgress(min=0, max=num_its)
display(f)

f1 = FloatProgress(min=0, max=k-1)
display(f1)

barpos = 1

for j in test_vals:
    results_array = np.array([])
    
    for i in range(0,k):
        t = completeTCPcalc(n = 100,
                        alphabeta_use = j[1],
                        alphabeta_sd_use = 2,
                        d = 2,
                        d_shift = j[0],
                        d_sd = 0,
                        n0 = 1E9,
                        max_d = 100,
                        dose_of_interest = j[2])
        results_array = np.append(results_array,TCP_cure_at_d_interest)
        f1.value = i    

        n = t[0]
        alphabeta_use = t[1]
        alphabeta_sd_use = t[2]
        d = t[3]
        d_shift = t[4]
        d_sd = t[5]
        n0 = t[6]
        TCP_pop = t[-2]
        TCPs = t[-3]
        nom_doses = t[-1]
        d_interest = t[8]
        frac_interest = t[9]
        TCP_cure_at_d_interest = t[10]
        max_d = 100

        results_array = np.append(results_array,TCP_cure_at_d_interest)

        ## include model parameters into array and insert before results
        param_array = np.array([])
        param_array = np.append(param_array,n)
        param_array = np.append(param_array,k)
        param_array = np.append(param_array,alphabeta_use)
        param_array = np.append(param_array,alphabeta_sd_use)
        param_array = np.append(param_array,d)
        param_array = np.append(param_array,d_shift)
        param_array = np.append(param_array,d_sd)
        param_array = np.append(param_array,n0)
        param_array = np.append(param_array,max_d)
        param_array = np.append(param_array,d_interest)
        #param_array = np.append(param_array,frac_of_interest)

    param_results_array = np.concatenate((param_array,results_array))
    saveasCSV(filename = 'TCPresults1.csv', array = param_results_array)

    #print(results_array)
    #print(param_array)
    
    barpos = barpos+1
    
    f.value = barpos # for updating progress bar
print(param_results_array[-1])


# In[ ]:



