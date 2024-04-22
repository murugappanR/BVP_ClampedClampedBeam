# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
from scipy.integrate import cumulative_trapezoid

E=169e9
t=0.5e-6
go=0.7e-6
d=1
l=80e-6
ep=8.854e-12
I=(1/12)*d*t**3
U=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
Y2=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
x=np.linspace(0,l,1400)

count=0
for u in U:    
   def bvp(x,W):
     w1, w2, w3, w4=W
     dw1=w2
     dw2=w3
     dw3=w4
     dw4=ep*d*U[count]**2/(2*E*I*(go-w1)**2)
     Q1=cumulative_trapezoid(w2**2,x,initial=0)
     Q2=(6/(l*t**2))*Q1[-1]*w3
     return [dw1,dw2,dw3,dw4]                          # add Q2 to dw4 to make it nonlinear
 
   def bc(Wa,Wb):
    w1a, w2a, w3a, w4a = Wa
    w1b, w2b, w3b, w4b = Wb
    return[w1a,w1b,
           w2a,w2b]

   a=1e-8;
   w1=a*np.sin(2*3.14*x/(2*l))
   w2=a*np.cos(2*3.14*x/(2*l))/(2*l)
   w3=-a*np.sin(2*3.14*x/(2*l))/(4*l**2)
   w4=-a*np.cos(2*3.14*x/(2*l))/(8*l**3)
   sol=solve_bvp(bvp, bc, x, [w1, w2, w3, w4])
   sol
   Y=sol.y[0]
   Y1=abs(Y)
   Y2[count]=Y1[700]
   plt.plot(x,sol.y[0])
   count=count+1 
#plt.plot(U,Y2)
