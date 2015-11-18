# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 12:58:36 2015

@author: Sara
"""
import numpy as np

GG_1 = np.array([[4,0,-4],
                 [0,4,-4],
                [4,4,8]])
BB_1 = np.array([[-5,0,5],
                 [0,-10,10],
                [-5,-10,-15]])
GG_2 = np.array([[4,0,-4],
                 [0,4,-4],
                [-4,-4,8]])
BB_2 = np.array([[-5,0,5], 
                 [0,-10,10],
                [5,10,-15]])

N = 3
vm = np.array([0,1.1249,0.93834])
an = np.array([0, 6.3, -3.44])
e = np.transpose(vm*np.cos(np.radians(an)))
f = np.transpose(vm*np.sin(np.radians(an)))

P_1 = [0]*N
Q_1 = [0]*N
P_2 = [0]*N
Q_2 = [0]*N
for i in range(0, N):
   s1_1 = 0.0
   s2_1 = 0.0 
   s1_2 = 0.0
   s2_2 = 0.0 
   for j in range(0, N):
       s1_1 = s1_1 + GG_1[i,j] * e[j] - BB_1[i,j] * f[j] 
       s2_1 = s2_1 + GG_1[i,j] * f[j] + BB_1[i,j] * e[j] 
       s1_2 = s1_2 + GG_2[i,j] * e[j] - BB_2[i,j] * f[j] 
       s2_2 = s2_2 + GG_2[i,j] * f[j] + BB_2[i,j] * e[j] 
   #print s1           
   P_1[i] = e[i] * s1_1 + f[i] * s2_1 
   Q_1[i] = f[i] * s1_1 - e[i] * s2_1
   P_2[i] = e[i] * s1_2 + f[i] * s2_2 
   Q_2[i] = f[i] * s1_2 - e[i] * s2_2


print P_2
print Q_2