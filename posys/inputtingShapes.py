# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 15:13:38 2016

@author: saf537
"""

import case14_mod
import pandas as pd
import pypower.api as pypo 
import numpy as np

load_shapes = pd.read_csv('loadShapes1.csv')
sample      = 100*load_shapes.iloc[0][load_shapes.columns[9:]]

ppc0      = case14_mod.case14_mod(busN = 1,dlt = 0, op_change=1) # trivial case: original solutin
ppopt0    = pypo.ppoption(PF_ALG=2,VERBOSE=0, OUT_ALL=0)
r0        = pypo.runpf(ppc0, ppopt0) 
n_trs     = sum(r0[0]['bus'][:,1] == 2)
n_bldgs   = sum(r0[0]['bus'][:,1] == 1)
transfs   = r0[0]['bus'][:,0][r0[0]['bus'][:,1]==2]
bdgs      = r0[0]['bus'][:,0][r0[0]['bus'][:,1]==1]
tm_max    = len(sample)
  # bus 8, power after solution

TS  = np.empty((n_bldgs,n_trs,tm_max))
initial_loads_vec = r0[0]['bus'][:,2]
cont = 0
for i in bdgs:
    tmp = 0
    ini_load = initial_loads_vec[i-1]
    for j in transfs:   
        for ld in sample:
            # NEED TO NORMALIZE CORRECTLY
            dw           = ini_load/ld - 1  
            ppc          = case14_mod.case14_mod(busN = i-1, dlt = dw, op_change=1)
            ppopt        = pypo.ppoption(PF_ALG=2, VERBOSE=0, OUT_ALL=0) 
            r            = pypo.runpf(ppc, ppopt)
            new_loads    = r[0]['gen'][:,2][1:]
            #ini_load     = r[0]['bus'][i-1][2]
            TS[cont,tmp] = new_loads[tmp]
            del(r)
        tmp += 1
    cont += 1            