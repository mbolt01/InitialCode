# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 16:54:36 2016

@author: mb22
"""

## NTCP calc start

## LKB Model

#%%
#import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import pandas as pd
import itertools

#%%
import TCP_NTCP
#from lmfit import Model, Parameter, fit_report

#%%
## calcaulte the cumulative sum of doses
doses = [2,2,1,2,2,2,2,2,1,2,3]

tot_dose = list(itertools.accumulate(doses))
print(doses)
print(tot_dose)
#%%

def td50_calc(td50_1,v,n):
    ## calcualte TD50(V) which takes into account the volume effect
    return td50_1/(v**n)

def u_calc(d,td50,m):
    ## calculation of the u value which is the limit for the NTCP integration
    return (d-td50)/(m*td50)
    
def ntcp_integrand(x):
    ## This is the part of the NTCP which is within the integral
    return sp.exp(-0.5*x**2)
    
def ntcp_calc(d,td50_1,v,m,n):
    """
    NTCP calculation based on LKB model.
    Parameters required are:
    d = total dose
    td50_1 = dose at which 50% of patients see effect for partial volume
    v = partial volume
    m = describes SD of TD50 SD~m.TD50(V). Inversly related to curve steepness.
    n = describes steepness of curve. n = 0 indicates a serial structure.
    """
    
    ## calcualte TD50(V)
    td50 = td50_calc(td50_1,v,n)
    #print(td50)
    
    ## calculate u for the integration limit
    u = u_calc(d,td50,m)
    #print(u)

    ## calculate NTCP value from input parameters
    ntcp = (1/(sp.sqrt(2*sp.pi)))*sp.integrate.quad(ntcp_integrand,-sp.inf,u)[0]
    #print(ntcp)
    
    return ntcp

def ntcp_fit_calc(dose_data, td50_1, v, m, n):
    ## calculate the NTCP at supplied dose points
    ntcp_fitted = []
    for dose in dose_data:
        ntcp_fitted.append(ntcp_calc(dose, td50_1, v, m, n))
    return ntcp_fitted
    
def ntcp_curve_calc(dose_range,td50_1, v, m, n):
    ## calcualte the ntcp curve up to a maximum dose
    ntcp_curve = []
    step = 0.1
    doses = np.arange(0,dose_range+step,step)
    for dose in doses:
        #doses.append(dose)
        ntcp_curve.append(ntcp_calc(dose, td50_1, v, m, n))
    return doses,ntcp_curve
    
##  fit the TD50_1,m and n parameters.

## dif in squares to minimise
def sum_square_difs(vals):
    ## must provide a list of tuples contianing mathcing pairs of (calc,ref)
    ## can combine lists like vals = list(zip(calcs,refs))
    ## this will be used in determining the differnce between fitted curve
    ## and data points
    
    return sum((vals[i][0] - vals[i][1])**2 for i in range(len(vals)))
    
## produce the effective volume calculation which might be useful.

##veff = sum(Di/D)^1/n*deltaVi

def veff_calc(cdvh_data,n):
    """
    Must provide:
    cdvh_data = cdvh as a list of tuples (dose,volume (%))
    n = describes parallel (n=1) or serial structures (n=0)
    
    Calculates on its way to veff:
    di = dose for each bin in differential dvh
    dmax = max dose to organ
    dvi = volume in specific dose bin i
    
    Returns:
    ddvh = tuple containing ddvh (dose, volume (%))
    veff = effective volume (%)
    """
    ## extract dose (di) and cdvh from input
    ## n cant be zero or calc fails due to divide by zero
    if n==0:
        n=0.0001
    di, cdvh = zip(*cdvh_data)
    
    ## calc ddvh
    ddvh = -np.gradient(cdvh) # note the negative sign
    
    ## find max dose = find where cdvh is first equal to 0
    dmax = di[cdvh.index(0)]

    ## dvi is equal to the ddvh at the point i. use to calc the bin values
    veff = sum([(ddvh[i]*((di[i]/dmax)**(1/n))) for i in range(len(di))])
   
    ## combine the dose and calced ddvh into single list of tuples
    ddvh_data = list(zip(di,ddvh))
    
    return veff,ddvh_data ## veff is a percentage value


def range_list(m, perc=None, dif=None,n=None, spacing=None):
    """
    A function to create a list of values of length 2n+1, or set spacing.
    n is the number of values either side of the mean to return
    The values are centred around the mean, m and 
    have a range extending from  +/- perc of m.
    values returned will not exceed the m+/-perc specified
    """
    ## ensure required parameters are passed
    
    if perc==None and dif==None:
        raise Exception('Need to specify a range with perc or dif')    
    if n==None and spacing==None:
        raise Exception('Need to specify number or spacing of output')
    if n!=None and spacing!=None:
        raise Exception('Ambiguous input as both n and spacing were supplied')
        
    ## convert percentage dif to absolute range
    if perc == None:
        abs_dif = dif
    if dif == None:
        abs_dif = m/100*perc
    #print(abs_dif)
            
    if spacing == None:
        if n < 1:
            if n == 0:
                results = [m]
            else:
                raise Exception('need at least 1 value either side of mean for n')
        else:
            n = np.floor(n) # wnat whole numbers either side
            results = np.linspace(m-abs_dif,m+abs_dif,2*n+1)
    
    if n==None:
        if spacing==0:
            results = [m]
        else:
            vals = []
            val = m
            ## add lower vlaues
            while val >= m-abs_dif:
                vals.append(val)
                val = val - spacing
                
            val = (m+spacing)
            ## add upper values
            while val <= m+abs_dif:
                vals.append(val)
                val = val + spacing
    
            results = sorted(vals)
    
    return list(results)

#%%

test = range_list(m=80,perc=10,spacing=2)
print(test)
#plt.plot(test,marker='o')
#plt.axhline(y=87.7)
#range_array(100,10,n=10)


#%%

## get doses for each fraction using fracdose function.
num_frac = 36
num_pts = 1000
for i in range(num_pts):
    alld = []
    for i in range(num_frac):
        a = TCP_NTCP.fracdose(2,0,0.5)
        alld.append(a)
    plt.plot(alld,alpha=0.05)

#%%

########### try and use lmfit to dot he fitting in a simpler way?
#
#ntcp_model = Model(ntcp_calc,independant_vars=['d'])
#
#model_result = ntcp_model.fit(ntcp_data,d=60,
#                              td50_1=Parameter(value=50, min=20, max=100),
#                              v=Parameter(value=1, vary=True, min=0, max=1),
#                              m=Parameter(value=0.5, min=0.1, max=1),
#                              n=Parameter(value=0.5, min=0.1, max=1))
#
#print(fit_report(model_result))
#
###########



def ntcp_data_fit(dose_data,ntcp_data,initial_params,v=1.0):
    """
    fucntion to fit teh NTCP model to supplied data and return the parameters.
    At somepoint in the process, if parameter values are not supplied
    this function will need calling to detemien them.
    i.e. if data is supplied, then fit the values, if not then use supplied vals.
    Funciton should only return fitted params, not do any plotting etc.
    """
    
    #plt.close() # close any open plots
    ## some example data to fit to and plot
    dose_data = dose_data#[55,60, 62, 67, 72, 65]
    ntcp_data = ntcp_data#[0.1,0.15,0.1,0.2,0.3, 0.19]
    
    ## specify some initial starting values
    initial_params = initial_params # supply inital params as a list to the function
    ## can supply all at once using *initial_params (must be in correct order)
        
    ## calculate NTCP for supplied data
    ntcp_fit = ntcp_fit_calc(dose_data,*initial_params)
    
    ## calc dif of squares (for use in optimisation)
    ntcp_dif_squares = sum_square_difs(list(zip(ntcp_data,ntcp_fit)))
    #print(ntcp_dif_squares)
    
    ## fit the parameters TD50_1, m, n using scipy
    ## note v_would be specified on a patient by patient basis in reality?
    ## but for my purposes could use fixed values to see the effect of changes?
    v_val_upper = v#1.0 # can set v in the function
    v_val_lower = v_val_upper*0.999
    
    set_bounds = ([0,v_val_lower,0,0],
                  [100,v_val_upper,1,1])
    
    #methods = ['dogbox','trf']
    ## could hold parameters fixed by specifying a very small range?
    
    all_results_list = []
    
    #for i in range(len(methods)):
        #print(methods[i])
    popt,pcov = sp.optimize.curve_fit(f = ntcp_fit_calc,
                                xdata = dose_data,
                                ydata = ntcp_data,
                                p0 = initial_params,
                                bounds = set_bounds,
                                method='trf') #method : {‘lm’, ‘trf’, ‘dogbox’}
    
    perr = np.sqrt(np.diag(pcov))
    #for i in range(len(popt)):
    #    labels = ['TD50_1: ', 'v: ','m: ', 'n: ']
    #    print(labels[i] + "%.2f" % popt[i])
    #print(results[0])
    #print(popt)
    #print(perr)
    
    ## calculate complete NTCP curve (using fitted params)
    #fitted_params = [param*1 for param in initial_params]
    fitted_params = [param for param in popt]
    fitted_params[1]=1
    max_dose = 100

    #doses, ntcp_curve = ntcp_curve_calc(max_dose,*fitted_params)
    
    #results_list = [methods[i],*fitted_params]
    #all_results_list.append(results_list)
        
        ## plot the results
        ## NTCP fitted points
        #plt.plot(dose_data,ntcp_fit,ls='',marker='^', label = 'initial guess points')
        ## NTCP curve
        #cols=['blue','red']
        #plt.plot(doses,ntcp_curve, label = methods[i], c=cols[i])
    
    ## plot the supplied data
    #plt.plot(dose_data,ntcp_data, ls='', marker='o', label = 'data')
    ## add a legend
    #plt.legend(loc='best', numpoints=1)
    #plt.show()
    #print('end')
    #print(all_results)
    
    #plt.close()
    #print(ntcp_curve[1])
    #abs_difs = [ntcp_curve[i]-ntcp_curve1[i] for i in range(len(ntcp_curve))]
    #ratio = [ntcp_curve[i]/ntcp_curve1[i] for i in range(len(ntcp_curve))]
    #print('dif')
    #plt.plot(ratio)
    #print(all_results_list)
    
    #columns = ['method', 'td50_1', 'v', 'm', 'n']
    #df = pd.DataFrame(all_results_list,columns=columns)
    #print(df)
    return popt # return the fitted params

    
#%%

## test of the ntcp function

## dose and ntcp data to use for fitting (lists)
d = [55,60, 62, 67, 72, 65]
nd = [0.1,0.15,0.1,0.3,0.3, 0.19]

## initial params (use these defaults, but allow list to be suppled)
intial_ntcp_td50_1 = 70
intial_ntcp_v = 1.0
intial_ntcp_m = 0.1
intial_ntcp_n = 0.1
#[td50,v,m,n)]
initial_params = [intial_ntcp_td50_1,intial_ntcp_v,intial_ntcp_m,intial_ntcp_n]

## optimise fitting parameters [td50,v,m,n] is returned
pop_fit = ntcp_data_fit(dose_data = d,
                     ntcp_data = nd,
                     initial_params = initial_params)

print(pop_fit)

## plot supplied data
plt.plot(d,nd, ls='', marker='o', label = 'data')

## plot fit - supply maximum dose and fitting parameters.
max_dose = 100
## could calculate range of each fitting parameter and provide to this to
## get multiple fits and then store these an an array/list
fit_x,fit_y = ntcp_curve_calc(max_dose,*pop_fit)
plt.plot(fit_x,fit_y,color='red')

## need to fit based on range of doses supplied






#%%



plt.close()
## some example data to fit to and plot
dose_data = [55,60, 62, 67, 72, 65]
ntcp_data = [0.1,0.15,0.1,0.2,0.3, 0.19]

## specify some initial starting values
intial_ntcp_td50_1 = 70
intial_ntcp_v = 1.0
intial_ntcp_m = 0.1
intial_ntcp_n = 0.1

initial_params = [intial_ntcp_td50_1,intial_ntcp_v,intial_ntcp_m,intial_ntcp_n]
## can supply all at once using *initial_params (must be in correct order)

## calc the ntcp at the supplied dose points
   
    
## calculate NTCP for supplied data
ntcp_fit = ntcp_fit_calc(dose_data,*initial_params)

## calc dif of squares (for use in optimisation)
ntcp_dif_squares = sum_square_difs(list(zip(ntcp_data,ntcp_fit)))
#print(ntcp_dif_squares)

## fit the parameters TD50_1, m, n using scipy
## note v_would be specified on a patient by patient basis in reality?
## but for my purposes could use fixed values to see the effect of changes?
v_val_upper = 1.0
v_val_lower = v_val_upper*0.999

set_bounds = ([0,v_val_lower,0,0],
              [100,v_val_upper,1,1])

methods = ['dogbox','trf']
## could hold parameters fixed by specifying a very small range?

all_results_list = []

for i in range(len(methods)):
    #print(methods[i])
    popt,pcov = sp.optimize.curve_fit(f = ntcp_fit_calc,
                                xdata = dose_data,
                                ydata = ntcp_data,
                                p0 = initial_params,
                                bounds = set_bounds,
                                method='trf') #method : {‘lm’, ‘trf’, ‘dogbox’}

    perr = np.sqrt(np.diag(pcov))
    #for i in range(len(popt)):
    #    labels = ['TD50_1: ', 'v: ','m: ', 'n: ']
    #    print(labels[i] + "%.2f" % popt[i])
    #print(results[0])
    #print(popt)
    #print(perr)

## calculate complete NTCP curve (using fitted params)
#fitted_params = [param*1 for param in initial_params]
    fitted_params = [param for param in popt]
    fitted_params[1]=1
    max_dose = 100

    doses, ntcp_curve = ntcp_curve_calc(max_dose,*fitted_params)
    
    results_list = [methods[i],*fitted_params]
    all_results_list.append(results_list)
    
    ## plot the results
    ## NTCP fitted points
    #plt.plot(dose_data,ntcp_fit,ls='',marker='^', label = 'initial guess points')
    ## NTCP curve
    cols=['blue','red']
    plt.plot(doses,ntcp_curve, label = methods[i], c=cols[i])

## plot the supplied data
plt.plot(dose_data,ntcp_data, ls='', marker='o', label = 'data')
## add a legend
plt.legend(loc='best', numpoints=1)
plt.show()
#print('end')
#print(all_results)

#plt.close()
#print(ntcp_curve[1])
#abs_difs = [ntcp_curve[i]-ntcp_curve1[i] for i in range(len(ntcp_curve))]
#ratio = [ntcp_curve[i]/ntcp_curve1[i] for i in range(len(ntcp_curve))]
#print('dif')
#plt.plot(ratio)
#print(all_results_list)

columns = ['method', 'td50_1', 'v', 'm', 'n']
df = pd.DataFrame(all_results_list,columns=columns)
print(df)

## the differnce between the absolute values for the 2 methods is small
## this is the case even if the change in paraemter values seems large.
## for my work i;m only interested in the shape of the curve, not the values.
## I can use actual values from the literature.

## When setting v ~ 1 trf does a better job of fitting the curve.
## note that when setting v=1 this effectively eliminates n as a parameter
## so only m controls the steepness.
## maybe this is the best method?
## both fitting algorithms give v sim results then.

#plt.show()

#%%

## do a comparison between the 2 fitting methods
a = df[df['method']=='dogbox']
b = df[df['method']=='trf']

list_results = df.values.tolist()
#print(df)

## get the numerical valeus vack from the df
a=list_results[0][1:]
b=list_results[1][1:]
print('a')
print(a)
print('b')
print(b)

difs = [a[i]-b[i] for i in range(len(a))]
print('dif')
print(difs)

## calc teh NTCP curve for these parameter sets

d1,ntcp1 = ntcp_curve_calc(100,*a)
d2,ntcp2 = ntcp_curve_calc(100,*b)
ntcp_dif = [((ntcp2[i]/ntcp1[i])-1)*100 for i in range(len(ntcp1))]

#plt.plot(d1,ntcp1)
#plt.plot(d2,ntcp2)
## another proof that differnt params can give almost identical curves

#plt.plot(d1,ntcp_dif)
#plt.ylabel('% dif between fits')
#plt.xlabel('Dose (Gy)')
## very small percentage difference between the fits

## how about the curve steepness?

#plt.plot(d1,ntcp1)
#plt.plot(d2,ntcp2)

d1_np = np.array(d1)
d2_np = np.array(d2)
ntcp1_np = np.array(ntcp1)
ntcp2_np = np.array(ntcp2)

#plt.plot(d1_np,ntcp1_np)
#plt.plot(d2_np,ntcp2_np)

ntcp1_grad = np.gradient(ntcp1)
ntcp2_grad = np.gradient(ntcp2)

plt.plot(d1,ntcp1_grad)
plt.plot(d2,ntcp2_grad)
plt.ylabel('gradient')
plt.xlabel('Dose (Gy)')
## also very minor differnces in the gradient

#%%
## demonstration oof different td50_1 and m params
## note v is held fixed at 1 so n plays no part

p,pp = ntcp_curve_calc(100,50,1,0.2,1)
q,qq = ntcp_curve_calc(100,55,1,0.3,0)

plt.plot(p,pp, label='p')
plt.plot(q,qq, label='q')
plt.legend(loc='best')

#%%

## might want to loop through all the pateints and calcualte veff for them?
## then will have a spread of the values to use.


## read in a file - use the large df I created? Then can loop through columns to get veff???
df_dvh = pd.read_excel(r'C:\Users\mb22\OneDrive\PhD\Quasar Shared\Data\Trials\PARSPORT\DVHs\To use\files-orig+sort+clean2+header_rename - Copy\1003.xls')

## take only l_parotid
df_l_par = df_dvh[['dose_(cgy)','l_parotid']]
## convert to Gy
df_l_par['dose_(cgy)']=df_l_par['dose_(cgy)']/100

## plot the cdvh
#plt.plot(df_l_par['dose_(cgy)'],df_l_par['l_parotid'])

## get dose and cdvh info
doses = df_l_par['dose_(cgy)'].values
cdvh = df_l_par['l_parotid'].values

## calc ddvh
ddvh = -np.gradient(cdvh)
## add back into df
df_l_par['ddvh']=ddvh

## create a list of tuples containing the cdvh data for use in function
dvh_data = list(zip(doses,cdvh))

dose_only, cdvh_only = zip(*dvh_data) # this unzips a zipped list

#print(dvh_data[:10])

#plt.plot(doses,ddvh)
#plt.plot(df_l_par['dose_(cgy)'],df_l_par['ddvh'])
#df_l_par.head()


#%%
#print(dvh_data)
evols = []
n_vals = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
for n in n_vals:
    eff_vol,ddvh_data = veff_calc(dvh_data,n)
    #plt.plot(ddvh)
    #print(ddvh_data)
    retn_dose, retn_ddvh = zip(*ddvh_data)
    #plt.plot(retn_dose, retn_ddvh)
    evols.append(eff_vol)
    print(eff_vol)
plt.plot(n_vals,evols)
plt.xlabel('n')
plt.ylabel('Effective Vol')
plt.xlim(0,1)
plt.show()
#plt.ylim(0.95,1.35)
#print(dvh_data)

#plt.plot(difs)

#plt.plot(dvh_data)