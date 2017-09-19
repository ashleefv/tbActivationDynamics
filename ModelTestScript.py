# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 23:02:44 2016

@author: Steve
"""


from RunModel import RunModel
import matplotlib.pyplot as plt
from WSDetModel import WSDetModel
from DetModel import DetModel
from ColDetModel import ColDetModel
import time

#Run model, and time how long it takes to calculate it
start = time.time()
modelToRun = ColDetModel()
y,t,p,sol,model = RunModel(model=modelToRun, inputMode='file', inputSource='ColDet-model-parameters.csv')
end = time.time()
print(end-start)

#Plot results
fontSize = 14
DPI = 160

plt.figure(1)
plt.semilogy(t, sol[:,15], 'c-', label = 'Intracellular Bacteria', linewidth = 3)
plt.semilogy(t, sol[:,16], 'k--', label = 'Extracellular Bacteria', linewidth = 3)
plt.semilogy(t, sol[:,19], 'g-.', label = 'Bacterial Leakage', linewidth = 3)
plt.axis([t[0], t[-1], 10**-2, 10**4])
plt.legend(loc=int(4), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Bacterial Count', fontsize = fontSize)
plt.savefig('bacteria_vs_time.jpg', dpi=DPI)
plt.show()
    
plt.figure(2)
plt.semilogy(t, sol[:,0], 'c-', label = 'Resting Macrophages', linewidth = 3)
plt.semilogy(t, sol[:,1] - sol[:,2], 'k--', label = 'Infected Macrophages', linewidth = 3)
plt.semilogy(t, sol[:,3], 'g-.', label = 'Activated Macrophages', linewidth = 3)
plt.axis([t[0],t[-1],1,10**6])
plt.legend(loc=int(7), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Macrophage Count', fontsize = fontSize)
plt.savefig('macrophages_vs_time.jpg', dpi=DPI)
plt.show()

plt.figure(3)
plt.plot(t, sol[:,4], 'c-', label = 'Th0 Cells', linewidth = 3)
plt.plot(t, sol[:,5], 'k--', label = 'Th1 Cells', linewidth = 3)
plt.plot(t, sol[:,6], 'g-.', label = 'Th2 Cells', linewidth = 3)
plt.plot(t, sol[:,7], 'b-', label = 'T80 Cells', linewidth = 3)
plt.plot(t, sol[:,8], 'r--', label = 'T8 Cells', linewidth = 3)
plt.plot(t, sol[:,9], 'm-.', label = 'Tc Cells', linewidth = 3)
#plt.legend(loc=int(4), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('T Cell Count', fontsize = fontSize)
plt.savefig('Tcells_vs_time.jpg', dpi=DPI)
plt.show()
    
plt.figure(4)
plt.plot(t, sol[:,10], 'c-', label = 'F_a', linewidth = 3)
plt.plot(t, sol[:,11], 'k--', label = 'I_y', linewidth = 3)
plt.plot(t, sol[:,12], 'g-.', label = 'I_12', linewidth = 3)
plt.plot(t, sol[:,13], 'b-', label = 'I_10', linewidth = 3)
plt.plot(t, sol[:,14], 'r--', label = 'I_4', linewidth = 3)
plt.legend(loc=int(7), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Cytokine Levels', fontsize = fontSize)
plt.savefig('cytokines_vs_time.jpg', dpi=DPI)
plt.show()

plt.figure(5)
plt.semilogy(t, sol[:,17], 'c-', label = 'Collagen', linewidth = 3)
plt.semilogy(t, sol[:,18], 'k--', label = 'MMP', linewidth = 3)
plt.legend(loc=int(4), fontsize = fontSize)
plt.axis([t[0],t[-1],10**-10,10**-2])
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Concentration (g$\cdot$cm$^{-3}$)', fontsize = fontSize)
plt.savefig('Collagen_and_MMP-1_vs_time.jpg', dpi=DPI)
plt.show()

plt.figure(6)
plt.plot(t, sol[:,19], 'c-', label = 'Bacterial Leakage', linewidth = 3)
plt.legend(loc=int(7), fontsize = fontSize)
#plt.axis([t[0],t[-1],1,10**4])
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Bacterial Count', fontsize = fontSize)
plt.savefig('Bacterial_Leakage.jpg', dpi=DPI)
plt.show()

plt.figure(7)
plt.plot(t, sol[:,10], 'c-', label = 'TNF', linewidth = 3)
plt.legend(loc=int(4), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Cytokine Levels', fontsize = fontSize)
plt.axis([t[0],t[-1],0,1])
plt.savefig('TNFalpha_vs_time.jpg', dpi=DPI)
plt.show()

