# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 10:37:31 2016

Using PYPOWER to estimate the impact of small variations in the 
obtained angles and voltages.

@author: saf537
"""



import pypower.api as pypo # to install, just pip install PYPOWER
import case14_mod



ppc   = case14_mod.case14_mod(busN = 4,dlt = 0.05)
ppopt = pypo.ppoption(PF_ALG=2) # Careful: have to use Newton Method!!!
r     = pypo.runpf(ppc, ppopt)
volt  = r[0]['bus'][r[0]['bus'][:,1] == 1,7]  # voltage 7, type 2 
ang   = r[0]['bus'][r[0]['bus'][:,1] != 3,8]  # angle 8, type 2
 
#pypo.printpf(r)