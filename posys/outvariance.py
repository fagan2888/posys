# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 10:37:31 2016

Using PYPOWER to estimate the impact of small variations in the 
obtained angles and voltages.
"""

import pypower.api as pypo # to install, just pip install PYPOWER
import case14_mod
from numpy import zeros, hstack, transpose
import matplotlib.pyplot as plt

#import matplotlib as matplot

dlt_vec   = [0.05, -0.05, 0.10, -0.10, 0.50, -0.50]
ppc0      = case14_mod.case14_mod(busN = 1,dlt = 0) # trivial case: original solutin
ppopt0    = pypo.ppoption(PF_ALG=2) # Careful: have to use Newton Method!!!
r0        = pypo.runpf(ppc0, ppopt0)       
m1        = r0[0]['bus'][r0[0]['bus'][:,1] == 1,7]
v0        = [m1 for i in range(0,8)] 
v0        = transpose(v0)
m2        = r0[0]['bus'][r0[0]['bus'][:,1] != 3,8] 
a0        = [m2 for i in range(0,9)] 
#a0        = transpose(a0)
mtr       = zeros([9,22,len(dlt_vec)]) 

j = 1
for j in range(0,len(dlt_vec)):
    for i in range(0,9):
        ppc        = case14_mod.case14_mod(busN = i,dlt = dlt_vec[j])
        ppopt      = pypo.ppoption(PF_ALG=2) # Careful: have to use Newton Method!!!
        r          = pypo.runpf(ppc, ppopt)
        volt       = r[0]['bus'][r[0]['bus'][:,1] == 1,7]  # voltage 7, type 2 
        ang        = r[0]['bus'][r[0]['bus'][:,1] != 3,8]  # angle 8, type 2
        mtr[i,:,j] = hstack((volt,ang))
        

for j in range(0,len(dlt_vec)):
    
    fig     = plt.figure()
    a       = fig.add_subplot(1,1,1)
    imgplot = plt.imshow(100*(mtr[:,:8,j]-v0) / v0,interpolation="nearest",
                         cmap="seismic")
    a.set_title('Voltages')
    vname   = 'Voltages' + str(dlt_vec[j]) + '.png'
    plt.colorbar(orientation ='horizontal')
    plt.xlabel('Number of loads')
    plt.ylabel('Load that changed')
    plt.savefig(vname, dpi=100)
    
    fig     = plt.figure()
    a=fig.add_subplot(1,1,1)
    imgplot = plt.imshow(100*(mtr[:,9:,j] - a0)/a0)
    a.set_title('Angles')
    aname   = 'Angles' + str(dlt_vec[j]) + '.png'
    plt.colorbar(orientation ='horizontal')
    plt.xlabel('Number of buses')
    plt.ylabel('Load that changed')
    plt.savefig(aname, dpi=100)

# 
#pypo.printpf(r)
    

