#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 13:54:27 2025

@author: chenhsiaoyi
"""

import matplotlib.pyplot as plt
import numpy as np

ang=[]
for i in range(0,200,20):
    ang.append(i)
#ang.append(89)

for i in ang:

    R=[]
    E=[]
    f = open("./R_tr_"+str(i), "r")
    #f = open("./Rij_trace_q1_rotM", "r")
    for line in f:
        a=line.split()
        if i ==0:
            E.append(float(a[0])-0.000)
        else:
            E.append(float(a[0]))
            
        R.append(float(a[2]))
    fig, ax = plt.subplots(figsize=(8, 1.3))
    
    plt.plot(E[:],R[:],color='blue')
    plt.xlim([-0.03,0.005])
    plt.ylim([-100.0,1600])
    #plt.xlim([-0.05,0.08])
    #plt.ylim([0.0,9])
    plt.savefig('R'+str(i)+'.pdf')
    
R=[]
E=[]

f = open("./R_tr_0", "r")
#f = open("./Rij_trace_q1_rotM", "r")
for line in f:
    a=line.split()
    E.append(float(a[0])*1000)
    R.append(float(a[2]))
fig, ax = plt.subplots(figsize=(8, 4))
plt.plot(E[:],R[:],color='blue')
plt.xlim([-30,5])
#plt.ylim([0.0,400])
#plt.xlim([-0.05,0.08])
plt.ylim([0.0,750])
plt.savefig('R0_large.pdf')
