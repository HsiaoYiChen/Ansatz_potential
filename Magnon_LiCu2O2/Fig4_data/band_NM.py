#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 13:20:51 2025

@author: chenhsiaoyi
"""

import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(6, 5))
   
mine_k=[]
mine_E=[]
k=[]
E=[]
n=0
f = open("./LiCu2O2_band.dat", "r")
for line in f:
    if (len(line) >3):
        a=line.split()
        #print(a,len(line) )
        k.append(float(a[0]))
        E.append(float(a[1]))
    else:
        n=n+1
        if (n >0) :
            mine_k.append(k)
            mine_E.append(E)
        k=[]
        E=[]
print(n)   
mine_k=np.asarray(mine_k)
mine_E=np.asarray(mine_E)
mine_E=mine_E-11.4
mine_k=mine_k/1.12
for i in range(len(mine_k[:,1])): 
    plt.plot(mine_k[i,:],mine_E[i,:],c='red',linestyle=':',zorder=10,lw=3)
    #plt.xlim([0,mine_k[-1,-1]])
    #plt.ylim([-1.5,1.5])

mine_k=[]
mine_E=[]
k=[]
E=[]
n=0
f = open("./LiCu2O2.dat.gnu_no_mag", "r")
for line in f:
    if (len(line) >2):
        a=line.split()
        #print(a,len(line) )
        k.append(float(a[0]))
        E.append(float(a[1]))
    else:
        n=n+1
        if (n >0) :
            mine_k.append(k)
            mine_E.append(E)
        k=[]
        E=[]
print(n)   
mine_k=np.asarray(mine_k)
mine_E=np.asarray(mine_E)
mine_E=mine_E-11.4
for i in range(len(mine_k[:,1])): 
    plt.plot(mine_k[i,:],mine_E[i,:],c='blue',zorder=3,lw=3)


plt.xlim([0,mine_k[-1,-1]])
plt.ylim([-1.5,1])
    
    
 
plt.xticks([])
#plt.hlines(0.0,0.0,mine_k[-1,-1],color='black',lw=1)
plt.vlines(mine_k[-1,20],-1.5,1.5,color='black',lw=1)
plt.vlines(mine_k[-1,40],-1.5,1.5,color='black',lw=1)
plt.vlines(mine_k[-1,60],-1.5,1.5,color='black',lw=1)
plt.vlines(mine_k[-1,80],-1.5,1.5,color='black',lw=1)
#plt.vlines(mine_k[-1,900],-1.5,1.5,color='black',lw=1)

#plt.vlines(mine_k[-1,1080],-1.5,1.5,color='black',lw=1)

#plt.vlines(mine_k[-1,1260],-1.5,1.5,color='black',lw=1)
#plt.vlines(mine_k[-1,1440],-1.5,1.5,color='black',lw=1)

plt.savefig('band_NM.pdf')
