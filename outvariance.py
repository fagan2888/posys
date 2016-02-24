# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 10:37:31 2016

Using PYPOWER to estimate the impact of small variations in the 
obtained angles and voltages.
"""

import pypower.api as pypo # to install, just pip install PYPOWER
import case14_mod
from numpy import zeros, hstack
import pandas as pd
from pylab import imshow
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
#import matplotlib as matplot

dlt_vec = [0.05, -0.05, 0.10, -0.10, 0.50, -0.50]
mtr     = zeros([9,22,len(dlt_vec)]) 

j = 1
for j in range(0,len(dlt_vec)):
    for i in range(0,9):
        ppc      = case14_mod.case14_mod(busN = i,dlt = dlt_vec[j])
        ppopt    = pypo.ppoption(PF_ALG=2) # Careful: have to use Newton Method!!!
        r        = pypo.runpf(ppc, ppopt)
        volt     = r[0]['bus'][r[0]['bus'][:,1] == 1,7]  # voltage 7, type 2 
        ang      = r[0]['bus'][r[0]['bus'][:,1] != 3,8]  # angle 8, type 2
        mtr[i,:,j] = hstack((volt,ang))
        

for j in range(0,len(dlt_vec)):
    
    fig = plt.figure()
    a=fig.add_subplot(1,1,1)
    imgplot = plt.imshow(mtr[:,:8,j])
    a.set_title('Voltages')
    vname = 'Voltages' + str(dlt_vec[j]) + '.png'
    plt.colorbar(orientation ='horizontal')
    plt.savefig(vname, dpi=100)
    
    fig = plt.figure()
    a=fig.add_subplot(1,1,1)
    imgplot = plt.imshow(mtr[:,9:,j])
    a.set_title('Angles')
    aname = 'Angles' + str(dlt_vec[j]) + '.png'
    plt.colorbar(orientation ='horizontal')
    plt.savefig(aname, dpi=100)

# 
#pypo.printpf(r)
    

