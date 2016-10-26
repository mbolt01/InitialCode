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
#from lmfit import Model, Parameter, fit_report

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

#%%
plt.close()
## some example data to fit to and plot
dose_data = [55,60, 62, 67, 72, 65]
ntcp_data = [0.1,0.15,0.1,0.2,0.3, 0.19]

## specify some initial starting values
intial_ntcp_td50_1 = 70
intial_ntcp_v = 1.0
intial_ntcp_m = 0.2
intial_ntcp_n = 0.2

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
                                method=methods[i]) #method : {‘lm’, ‘trf’, ‘dogbox’}

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