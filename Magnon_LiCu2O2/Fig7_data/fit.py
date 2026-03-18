#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 21:48:13 2025

@author: chenhsiaoyi
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


fit_number=4
k=[]
E=[]
i=0
f = open("data", "r")
E_tmp=[]
for line in f:
    a=line.split()
    if len(a)>0:
        k.append(float(a[1]))
        t=[float(a[2]),float(a[3]),float(a[4]),float(a[5]),float(a[6]),float(a[7])]
        E.append(t)


k=np.asarray(k)
E=np.asarray(E)
E0=E[0,:]
E0=np.array(E0)*0.0

E[:,0]=E[:,0]-E[0,0]
E[:,1]=E[:,1]-E[0,1]
E[:,2]=E[:,2]-E[0,2]
E[:,3]=E[:,3]-E[0,3]
E[:,4]=E[:,4]-E[0,4]
E[:,5]=E[:,5]-E[0,5]
color=['blue','green','grey','orange','red','purple']
fig, ax = plt.subplots(figsize=(3, 8))
# Define the model function
def model_func(x, a, b,c,d):
    return a * x**2 + b * x**4+c* x**6+ d * x**8

# Fit the model to the data
for fit_number in [0,1,2,3,4,5]:

    params, covariance = curve_fit(model_func, k, E[:,fit_number])
    
# Extract fitted parameters
    a_fit, b_fit,c_fit,d_fit = params
    print(f"Fitted parameters: a = {a_fit}, b = {b_fit}")

# Optional: plot data and fitted curve
    x_fit = np.linspace(min(k), max(k), 500)
    y_fit = model_func(x_fit, a_fit, b_fit,c_fit,d_fit)


    plt.scatter(k, E[:,fit_number]+E0[fit_number], label="Data",c=color[fit_number])
    plt.plot(x_fit, y_fit+E0[fit_number], color=color[fit_number], label="Fit: $ax^2+bx^4+cx^6+dx^8$")
    min_ind=np.argmin(y_fit)
    print(x_fit[int(min_ind)])
#plt.legend()
#plt.xlabel("qy")
#plt.ylabel("energy(eV)")
#plt.title("Curve Fit")
plt.xlim([0.0,0.375])

plt.savefig('band_fit.pdf')
