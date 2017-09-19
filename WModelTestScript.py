# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 21:03:06 2016

@author: Steve
"""

from RunModel import RunModel
import matplotlib.pyplot as plt
from WDetModel import WDetModel

y,t,p,sol = RunModel(ModelName=WDetModel, inputMode='file', inputSource='Wigginton-model-parameters.csv')

plt.figure(1)
#plt.semilogy(t, sol[:,10], 'c-', label = 'B_I')
#plt.semilogy(t, sol[:,11], 'k--', label = 'B_E')
plt.plot(t, sol[:,10]+sol[:,11], 'g-.', label = 'Total bacteria')
plt.legend(loc=int(7))
plt.show()
    
plt.figure(2)
#plt.semilogy(t, sol[:,0], 'c-', label = 'M_R')
#plt.semilogy(t, sol[:,1], 'k--', label = 'M_I')
plt.plot(t, sol[:,2], 'g-.', label = 'M_A')
#plt.axis([t[0],t[-1],1,10**6])
plt.legend(loc=int(7))
plt.show()

plt.figure(3)
plt.plot(t, sol[:,3], 'c-', label = 'T_0')
plt.plot(t, sol[:,4], 'k--', label = 'T_1')
plt.plot(t, sol[:,5], 'g-.', label = 'T_2')
plt.plot(t, sol[:,3]+sol[:,4]+sol[:,5], 'r-.', label = 'total T cells')

plt.legend(loc=int(7))
plt.show()
    
plt.figure(4)
plt.plot(t, sol[:,6], 'k--', label = 'I_y')
plt.plot(t, sol[:,7], 'g-.', label = 'I_12')
plt.plot(t, sol[:,8], 'b-', label = 'I_10')
plt.plot(t, sol[:,9], 'r--', label = 'I_4')
plt.legend(loc=int(7))
plt.show