#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 13:54:27 2025

@author: chenhsiaoyi
"""

import matplotlib.pyplot as plt
import numpy as np

q_list=[91,85,79,73,67,61,55,49,43,37,31,25,19,13,7,1,193,
        385,577,769,775,781,787,793,799,805,811,817,823,829,835,841,847,853,859,865]

fig, ax = plt.subplots(figsize=(4.5, 5))
   # 1, 1, figsize=(8, 12), subplot_kw={'projection': '3d'})
max_R=0.0
min_R=0.0
for i in range(len(q_list)):

    R=[]
    E=[]
    f = open("./Rij_trace_q"+str(q_list[i])+"_12", "r")
    #f = open("./Rij_trace_q1_rotM", "r")
    count=1
    for line in f:
        a=line.split()
        if np.remainder(count,1) ==0:
            if (max_R < np.log(float(a[2]))):
                max_R =np.log(float(a[2]))
            if (min_R > np.log(float(a[2]))):
                min_R =np.log(float(a[2]))    
            E.append(float(a[0]))    
            R.append(np.log(float(a[2])))
        count=count+1
    print(min_R,max_R)
    q=[i]*len(R)
    R[0]=7.7528
    R[1]=-4.27970357541
    plt.scatter(q,E,c=R,marker='s')
    
    
    #plt.xlim([-0.02,0.01])
    plt.ylim([-0.04,0.08])
    #plt.xlim([-0.05,0.08])
    #plt.ylim([0.0,9])
plt.colorbar()
plt.savefig('R1.jpg')



