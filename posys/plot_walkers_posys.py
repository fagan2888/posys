# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 12:18:42 2016

@author: UO 

"""
import matplotlib.pyplot as plt
import corner
import numpy as np

median_20_init0,median_20_init1,median_20_init2,median_20_init3,median_20_init4,mediacn_50_init0,median_50_init1,median_50_init2,median_50_init3,median_50_init4 = np.load('random/median_300Walk9build500steps20.0valRandomInit0.npy'),np.load('random/median_300Walk9build500steps20.0valRandomInit1.npy'),np.load('random/median_300Walk9build500steps20.0valRandomInit2.npy'),np.load('random/median_300Walk9build500steps20.0valRandomInit3.npy'),np.load('random/median_300Walk9build500steps20.0valRandomInit4.npy'),np.load('random/median_300Walk9build500steps50.0valRandomInit0.npy'), np.load('random/median_300Walk9build500steps50.0valRandomInit1.npy'), np.load('random/median_300Walk9build500steps50.0valRandomInit2.npy'),np.load('random/median_300Walk9build500steps50.0valRandomInit3.npy'),np.load('random/median_300Walk9build500steps50.0valRandomInit4.npy')

m_inits = np.median(np.median([median_20_init0,median_20_init1,median_20_init2,median_20_init3,median_20_init4,mediacn_50_init0,median_50_init1,median_50_init2,median_50_init3,median_50_init4],1),1)
m_bdgs = np.median(np.median([median_20_init0,median_20_init1,median_20_init2,median_20_init3,median_20_init4,mediacn_50_init0,median_50_init1,median_50_init2,median_50_init3,median_50_init4],1),0)


fl_20_init0,fl_20_init1,fl_20_init2,fl_20_init3,fl_20_init4,fl_50_init0,fl_50_init1,fl_50_init2,fl_50_init3,fl_50_init4 = np.load('random/300Walk9build500steps20.0valRandomInit0.npy'),np.load('random/300Walk9build500steps20.0valRandomInit1.npy'),np.load('random/300Walk9build500steps20.0valRandomInit2.npy'),np.load('random/300Walk9build500steps20.0valRandomInit3.npy'),np.load('random/300Walk9build500steps20.0valRandomInit4.npy'),np.load('random/300Walk9build500steps50.0valRandomInit0.npy'), np.load('random/300Walk9build500steps50.0valRandomInit1.npy'), np.load('random/300Walk9build500steps50.0valRandomInit2.npy'),np.load('random/300Walk9build500steps50.0valRandomInit3.npy'),np.load('random/300Walk9build500steps50.0valRandomInit4.npy')

# -- Plotting example

plt.plot((fl_50_init1[:,:,2]).T)

"""
ndim, nwalkers = 9, 300
nsteps = 500
v = [1,2,3,4,5]

theta  = [v[0]*np.ones(ndim)]#ppc0["bus"][ind,2].copy()#+ 1e-0*np.random.randn(ndim)
theta *= theta>0.0

Walk9build500steps1,Walk9build500steps2,Walk9build500steps5, Walk9build500steps10 ,Walk9build500steps50 = np.load('300Walk9build500steps1.0val.npy'),np.load('300Walk9build500steps2.0val.npy'),np.load('300Walk9build500steps5.0val.npy'),np.load('300Walk9build500steps10.0val.npy'), np.load('300Walk9build500steps50.0val.npy') 
vardic = {Walk9build500steps1:1, Walk9build500steps2:2, Walk9build500steps5:5, Walk9build500steps10:10 , Walk9build500steps50:50}
varlist = [Walk9build500steps1,Walk9build500steps2,Walk9build500steps5, Walk9build500steps10 ,Walk9build500steps50]

samples = Walk9build500steps1[:, :, :].reshape((-1, ndim))

f, ((ax1, ax2, ax3, ax4, ax5)) = plt.subplots(len(v), ndim, sharex='col', sharey='row')
dic = {str(var): a for var in vardic.keys() for a in (ax1, ax2, ax3, ax4, ax5)}
for va in varlist:
    for nd in range(1,ndim+1):
        dic[int(va)].plot(va[nd,:,:])

fig = corner.corner(samples, labels=["$b1$", "$b2$", "$b3$","$b4$","$b5$",
                                     "$b6$","$b7$","$b8$","$b9$"],
                    truths=theta)
"""
