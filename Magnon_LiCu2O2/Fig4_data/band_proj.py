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
mine_w=[]
k=[]
E=[]
w=[]
n=0
f = open("./LiCu2O2_band.dat_dxy", "r")
for line in f:
    if (len(line) >4):
        a=line.split()
        #print(a,len(line) )
        k.append(float(a[0]))
        E.append(float(a[1]))
        w.append(float(a[2]))
    else:
        n=n+1
        if (n >50) :
            mine_k.append(k)
            mine_E.append(E)
            mine_w.append(w)
        k=[]
        E=[]
        w=[]
print(n)   
mine_k=np.asarray(mine_k)
mine_E=np.asarray(mine_E)
mine_w=np.asarray(mine_w)
mine_E=mine_E-11.4
mine_k=mine_k/1.12

mine_w[0,0]=1
plt.scatter(mine_k[:,:],mine_E[:,:],c=mine_w[:,:],linestyle=':',zorder=10,lw=3)
#plt.colorbar()
#for i in range(len(mine_k[:,1])): 
#    plt.scatter(mine_k[i,:],mine_E[i,:],c=mine_w[i,:],linestyle=':',zorder=10,lw=3)

plt.xlim([0.0,2.85])
plt.ylim([-1.5,0.5])
plt.savefig('proj_dxy.jpg')

fig, ax = plt.subplots(figsize=(6, 5))


mine_k=[]
mine_E=[]
mine_w=[]
k=[]
E=[]
w=[]
n=0
f = open("./LiCu2O2_band.dat_px", "r")
for line in f:
    if (len(line) >4):
        a=line.split()
        #print(a,len(line) )
        k.append(float(a[0]))
        E.append(float(a[1]))
        w.append(float(a[2]))
    else:
        n=n+1
        if (n >50) :
            mine_k.append(k)
            mine_E.append(E)
            mine_w.append(w)
        k=[]
        E=[]
        w=[]
print(n)   
mine_k=np.asarray(mine_k)
mine_E=np.asarray(mine_E)
mine_w=np.asarray(mine_w)
mine_E=mine_E-11.4
mine_k=mine_k/1.12

mine_w[0,0]=1
plt.scatter(mine_k[:,:],mine_E[:,:],c=mine_w[:,:],linestyle=':',zorder=10,lw=3)
#plt.colorbar()
plt.xlim([0.0,2.85])
plt.ylim([-1.5,0.5])
plt.savefig('proj_px.jpg')




fig, ax = plt.subplots(figsize=(6, 5))


mine_k=[]
mine_E=[]
mine_w=[]
k=[]
E=[]
w=[]
n=0
f = open("./LiCu2O2_band.dat_py", "r")
for line in f:
    if (len(line) >4):
        a=line.split()
        #print(a,len(line) )
        k.append(float(a[0]))
        E.append(float(a[1]))
        w.append(float(a[2]))
    else:
        n=n+1
        if (n >50) :
            mine_k.append(k)
            mine_E.append(E)
            mine_w.append(w)
        k=[]
        E=[]
        w=[]
print(n)   
mine_k=np.asarray(mine_k)
mine_E=np.asarray(mine_E)
mine_w=np.asarray(mine_w)
mine_E=mine_E-11.4
mine_k=mine_k/1.12

mine_w[0,0]=1
plt.scatter(mine_k[:,:],mine_E[:,:],c=mine_w[:,:],linestyle=':',zorder=10,lw=3)
#plt.colorbar()

plt.xlim([0.0,2.85])
plt.ylim([-1.5,0.5])
plt.savefig('proj_py.jpg')




