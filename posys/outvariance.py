# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 10:37:31 2016

Using PYPOWER to estimate the impact of small variations in the 
obtained angles and voltages.
"""

import os
import pypower.api as pypo # to install, just pip install PYPOWER
import case14_mod
from numpy import zeros, hstack, transpose
import matplotlib.pyplot as plt

#import matplotlib as matplot

dlt_vec   = [0.05, -0.05, 0.10, -0.10, 0.50, -0.50]
ppc0      = case14_mod.case14_mod(busN = 1,dlt = 0) # trivial case: original solutin
#ppopt0    = pypo.ppoption(PF_ALG=2, VERBOSE=0,
#                          OUT_ALL      = False,
#                          OUT_ALL_LIM  = False,
#                          OUT_AREA_SUM = False,
#                          OUT_BRANCH   = False,
#                          OUT_BUS      = False,
#                          OUT_GEN      = False,
#                          OUT_LINE_LIM = False,
#                          OUT_PG_LIM   = False,
#                          OUT_QG_LIM   = False,
#                          OUT_SYS_SUM  = False,
#                          OUT_V_LIM    = False) # Careful: have to use Newton Method!!!
ppopt0    = pypo.ppoption(PF_ALG=2, VERBOSE=0, OUT_ALL=0) # Careful: have to use Newton Method!!!
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
        ppopt      = pypo.ppoption(PF_ALG=2, VERBOSE=0, OUT_ALL=0) # Careful: have to use Newton Method!!!
        r          = pypo.runpf(ppc, ppopt)
        volt       = r[0]['bus'][r[0]['bus'][:,1] == 1,7]  # voltage 7, type 2 
        ang        = r[0]['bus'][r[0]['bus'][:,1] != 3,8]  # angle 8, type 2
        mtr[i,:,j] = hstack((volt,ang))
        

# -- plot results
figs = []
axs = []
ims = []

for j in range(0,len(dlt_vec)):

    figs.append(plt.figure())
    axs.append(figs[2*j].add_subplot(1,1,1))
    ims.append(axs[2*j].imshow(100*(mtr[:,:8,j]-v0) / v0,
                              interpolation="nearest",
                               cmap="seismic"))
    axs[2*j].set_title('Voltages (Power change: {0})'.format(dlt_vec[j]))
    vname   = 'Voltages' + str(dlt_vec[j]) + '.png'
    figs[2*j].colorbar(ims[2*j],orientation ='horizontal')
    axs[2*j].set_xlabel('Number of loads')
    axs[2*j].set_ylabel('Load that changed')
    figs[2*j].savefig(os.path.join("../output",vname), dpi=100, clobber=True)
    
    figs.append(plt.figure())
    axs.append(figs[2*j+1].add_subplot(1,1,1))
    ims.append(axs[2*j+1].imshow(100*(mtr[:,9:,j]-a0) / a0,
                              interpolation="nearest",
                               cmap="seismic"))
    axs[2*j+1].set_title('Angles (Power change: {0})'.format(dlt_vec[j]))
    aname   = 'Angles' + str(dlt_vec[j]) + '.png'
    figs[2*j+1].colorbar(ims[2*j+1],orientation ='horizontal')
    axs[2*j+1].set_xlabel('Number of buses')
    axs[2*j+1].set_ylabel('Load that changed')
    figs[2*j+1].savefig(os.path.join("../output",aname), dpi=100, clobber=True)
