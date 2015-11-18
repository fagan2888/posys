# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 20:24:20 2015

@author: Sara
"""
import numpy as np



jimag = np.sqrt(-1)
nline = len(ns)
nbus = len()
y   =  1 / z 

% initialize an empty sparse Ybus
Y = sparse(nbus,nbus) ; 

% buid Y bus
for i = 1:nline
  if tr(i) == 0 ;     % this line has no transformer
     tr(i) =  1 ;
  end
  tps      = tr(i)*exp(jimag*ts(i)*pi/180);
  from     = ns(i) ;
  to       = nr(i) ;
  Y(from,to  ) = Y(from,to  ) -  y(i) / conj(tps);
  Y(to  ,from) = Y(to  ,from) -  y(i) / tps;
  Y(from,from) = Y(from,from) + (y(i) + jimag*bc(i)/2) / (tps*conj(tps));
  Y(to  ,  to) = Y(to  ,  to) +  y(i) + jimag*bc(i)/2;
end 

% include the bus to ground conductance and susceptance
Gb = bus(:,8);         % bus conductance
Bb = bus(:,9);         % bus susceptance
if ((Gb ~=0) | (Bb ~=0) )
  Yb = diag(Gb+jimag*Bb);     
   Y = Y + Yb;
end

% separate real and imaginary parts since my power flow 
% programs compute the Jacobian and powers using them

GG = real(Y) ;
BB = imag(Y) ;

return
