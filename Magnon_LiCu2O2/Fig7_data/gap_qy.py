#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 22:19:15 2025

@author: chenhsiaoyi
"""

import numpy as np
import matplotlib.pyplot as plt
q=[0.1202,0.1593,0.1788,0.1969,0.2134,0.21943]
q2=[0.25]*6
scale=[2.15,1.8,1.6,1.4,1.2,1.0]
q=np.asarray(q[:])
q=1.0/q
q2=np.asarray(q2[:])
q2=1.0/q2
gap=[0.6532,0.449,0.3504,0.2387,0.1286,-0.010]
gap=np.asarray(gap[:])
gap[:]=gap[:]*1000

extra_q=[1/0.1202,12.5]
extra_g=[653.2,1000.0]
extra_q2=[4.0,4.0]

fig, ax = plt.subplots(figsize=(8, 5))

#plt.hlines(6, 0,1000)
#plt.vlines(400, 0,12)
plt.plot(gap,q,c='#1f77b4',lw=4)
#plt.plot(gap,q2,c='#ff7f0e',lw=4)
plt.plot(extra_g,extra_q,linestyle='--',c='#1f77b4',lw=2)
#plt.plot(extra_g,extra_q2,linestyle='--',c='#ff7f0e',lw=2)
#plt.xlim([-280,1050])
plt.ylim([3.9,12.7])
plt.savefig('gap_qy.pdf')
plt.show()


width=[0.06632,0.0715, 0.07513	,0.07997	,0.08759,0.09589		]
width=np.asarray(width)*1000

width2=[0.02674,0.03402, 0.0397	,0.04676	,0.05636,0.06627		]
width2=np.asarray(width2)*1000

extra_width=[0.06632*1000,0.05832*1000]
extra_width2=[0.02674*1000,0.01674*1000]

fig, ax = plt.subplots(figsize=(8, 3))
plt.plot(gap,width,c='#1f77b4',lw=2)
plt.plot(extra_g,extra_width,linestyle='--',c='#1f77b4')
plt.plot(gap,width2,c='#ff7f0e',lw=2)
plt.plot(extra_g,extra_width2,linestyle='--',c='#ff7f0e')
plt.savefig('gap_width.pdf')