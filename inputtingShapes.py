# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 15:13:38 2016

@author: saf537
"""

import case14_mod
import pandas as pd
import pypower.api as pypo 

load_shapes = pd.read_csv('loadShapes1.csv')
sample      = 10*load_shapes.iloc[0][load_shapes.columns[9:]]

ppc0      = case14_mod.case14_mod(busN = 1,dlt = 0) # trivial case: original solutin
ppopt0    = pypo.ppoption(PF_ALG=2) # Careful: have to use Newton Method!!!
r0        = pypo.runpf(ppc0, ppopt0) 
  # bus 8, power after solution

for i in range(0,14):
    ini_load = r0[0]['bus'][i][2]
    for ld in sample:
        dw  = ld / ini_load - 1    
        ppc = case14_mod.case14_mod(busN = i,dlt = dw)
        ppopt     = pypo.ppoption(PF_ALG=2) # Careful: have to use Newton Method!!!
        r         = pypo.runpf(ppc0, ppopt0)
        new_loads = r[0]['bus'][:,2]
        
        