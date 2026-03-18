#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 10:47:23 2025

@author: chenhsiaoyi
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

red_n=2
k_list=['0.0000','0.1667','0.2000','0.3333','0.4000','0.5000']
data=[]
weight=[1.0/5+1.0/6,1.0/6,1.0/5,1.0/6,1.0/5,1.0/6]

for k in k_list:
    f = open("R_ij_k_"+k, "r")
    i=red_n
    for line in f:
        a=line.split()
        if i ==red_n:
            data.append([float(k),float(a[0]),float(a[2])])
            i=0
        i=i+1
        
data=np.array(data)
fig, ax = plt.subplots(figsize=(1.3, 5.5))
plt.scatter(data[:,0],data[:,1],c=np.log(data[:,2]),marker='s',s=12)
ax.set_facecolor('gainsboro')

exp_data=[]
f = open("exp_data", "r")
for line in f:
    a=line.split()
    exp_data.append([float(a[0]),float(a[1])])
exp_data=np.array(exp_data)
exp_data[:,0]=exp_data[:,0]*6
exp_data[:,1]=exp_data[:,1]/1000
#plt.plot(exp_data[:,0]-exp_data[0,0],exp_data[:,1],c='red')
#plt.plot(exp_data[:,0]-exp_data[0,0]-0.5,exp_data[:,1],c='red')
plt.ylim([-0.0001,0.06])
plt.xlim([-0.01,0.51])
plt.savefig('Sqw.jpg')
plt.show()


k_list=['0.0000','0.1667','0.2000','0.3333','0.4000','0.5000']
data=[]
Aw=[0.0]*320
Aw=np.array(Aw)
j=0
w=[]
for k in k_list:
    f = open("R_ij_k_"+k, "r")
    i=0
    for line in f:
        a=line.split()
        if j==0:
            w.append(float(a[0]))
        Aw[i]=Aw[i]+float(a[2])*weight[j]
        i=i+1
    j=j+1
        






plt.plot(w,Aw)
plt.xlim([-0.0001,0.06])
plt.savefig('Sqw_sum.jpg')

fig, ax = plt.subplots(figsize=(6, 8))
simple_data=[]
f = open("simple", "r")
for line in f:
    a=line.split()
    simple_data.append([float(a[0]),float(a[1])])
simple_data=np.array(simple_data)
ax.fill_between(exp_data[:,0]-exp_data[0,0],exp_data[:,1]-0.0005,exp_data[:,1]+0.0005,color='pink', alpha=0.2)
plt.plot(exp_data[:,0]-exp_data[0,0],exp_data[:,1],lw=3,c='red')
plt.plot(exp_data[:,0]-exp_data[0,0],exp_data[:,1]+0.0005,c='pink')
plt.plot(exp_data[:,0]-exp_data[0,0],exp_data[:,1]-0.0005,c='pink')
plt.plot(-(exp_data[:,0]-exp_data[0,0])+1.0,exp_data[:,1],lw=3,c='red')
plt.plot(-(exp_data[:,0]-exp_data[0,0])+1.0,exp_data[:,1]+0.0005,c='pink')
plt.plot(-(exp_data[:,0]-exp_data[0,0])+1.0,exp_data[:,1]-0.0005,c='pink')
ax.fill_between(-(exp_data[:,0]-exp_data[0,0])+1.0,exp_data[:,1]-0.0005,exp_data[:,1]+0.0005,color='pink', alpha=0.2)

#plt.plot(exp_data[:,0]-exp_data[0,0]-2.07,exp_data[:,1],c='red')
#plt.plot(exp_data[:,0]-exp_data[0,0]-2.07,exp_data[:,1]+0.0005,c='grey')
#plt.plot(exp_data[:,0]-exp_data[0,0]-2.07,exp_data[:,1]-0.0005,c='grey')
#ax.fill_between(exp_data[:,0]-exp_data[0,0]-2.07,exp_data[:,1]-0.0005,exp_data[:,1]+0.0005,color='grey', alpha=0.2)





k_data=[0.1667,0.2000,0.3333,0.4000,0.5000]
E_data=[0.000659,0.001677,0.0034495,0.00447,0.00630]
initial_guess = [0.01, 0.01]
k_data=np.asarray(k_data)
E_data=np.asarray(E_data)
def F_model(q, A, B):
    return A + B*(q-1.0)**2

popt, pcov = curve_fit(F_model, k_data, E_data, p0=initial_guess)
A_fit, B_fit = popt
print("Fitted parameters:")
print(f"A = {A_fit:.4f}, B = {B_fit:.4f}")
perr = np.sqrt(np.diag(pcov))
print(f"Uncertainty in A = {perr[0]:.4f}")
print(f"Uncertainty in B = {perr[1]:.4f}")
q_fit = np.linspace(0.0, 1.0, 500)
F_fit = F_model(q_fit, *popt)

popt_low=A_fit-0.0004, B_fit-0.0008
F_fit_low = F_model(q_fit, *popt_low)

popt_high=A_fit+0.0004, B_fit+0.0008
F_fit_high = F_model(q_fit, *popt_high)

plt.plot(q_fit, F_fit, 'r-', lw=3, label="Fitted function",color='deepskyblue',zorder=1)
plt.plot(q_fit, F_fit_high, 'r-', lw=2, label="Fitted function",color='cyan',zorder=1)
plt.plot(q_fit, F_fit_low, 'r-', lw=2, label="Fitted function",color='cyan',zorder=1)
ax.fill_between(q_fit,F_fit_low,F_fit_high,color='cyan', alpha=0.2)


popt_low=A_fit-0.0004, B_fit-0.0008
F_fit_low = F_model(q_fit, *popt_low)


plt.plot(-q_fit+1.0, F_fit, 'r-', lw=3, label="Fitted function",color='deepskyblue',zorder=1)
plt.plot(-q_fit+1.0, F_fit_high, 'r-', lw=2, label="Fitted function",color='cyan',zorder=1)
plt.plot(-q_fit+1.0, F_fit_low, 'r-', lw=2, label="Fitted function",color='cyan',zorder=1)
ax.fill_between(-q_fit+1.0,F_fit_low,F_fit_high,color='cyan', alpha=0.2)




plt.scatter(simple_data[:,0],simple_data[:,1],c='blue',s=90,zorder=10)
plt.ylim([-0.0001,0.010])
plt.xlim([-0.001,0.5])
plt.savefig('Sqw_simple.pdf')
plt.show()
