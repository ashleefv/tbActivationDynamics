# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 15:38:18 2017

@author: Steve
"""


from RunModel import RunModel
import matplotlib.pyplot as plt
from WSDetModel import WSDetModel
from DetModel import DetModel
from ColDetModel import ColDetModel
import time
import numpy as np
import types


modelToRun = ColDetModel()
y,t,p,sol,model = RunModel(model=modelToRun, inputMode='file', inputSource='ColDet-model-parameters.csv')

TcellDepletion = [sol[:][-1], np.arange(1500,10001),p]
#TcellDepletion[0][4] = 0
TcellDepletion[2]['a_1a'] = 0
model=ColDetModel()

model.baseParameters = dict(p)
model.parameterAdjustFuncs['a_2'] = lambda y,t,p,pb: np.max([pb['a_2'] * np.exp(-(t-1500)/3000), pb['a_2']* np.exp(-(10000-1500)/3000)])
model.parameterAdjustFuncs['sr_1b'] = lambda y,t,p,pb: np.max([pb['sr_1b'] * np.exp(-(t-1500)/4500), pb['sr_1b']* np.exp(-(10000-1500)/4500)])

#model.parameterAdjustFuncs['a_2'] = lambda y,t,p,pb: pb['a_2'] * np.exp(-(t-1500)/500)
#model.parameterAdjustFuncs['sr_1b'] = lambda y,t,p,pb: pb['sr_1b'] * np.exp(-(t-1500)/750)
#model.parameterAdjustFuncs['sr_3b'] = lambda y,t,p,pb: pb['sr_3b'] * np.exp(-(t-1500)/750)
#model.parameterAdjustFuncs['sr_3b2'] = lambda y,t,p,pb: pb['sr_3b2'] * np.exp(-(t-1500)/750)
#model.parameterAdjustFuncs['sr_3bc'] = lambda y,t,p,pb: pb['sr_3bc'] * np.exp(-(t-1500)/750)


#def dT_0dt(self,y,t,p):
#         
#        return 0
#        
#model.dT_0dt = types.MethodType(dT_0dt, model)


y2,t2,p2,sol2,model2 = RunModel(model=model, inputMode='arg', paramsIn = TcellDepletion)

t = np.append(t,t2)
sol = np.append(sol, sol2, 0)

fontSize = 14
DPI = 160

plt.figure(1)
plt.semilogy(t, sol[:,15], 'c-', label = 'Intracellular Bacteria', linewidth = 3)
plt.semilogy(t, sol[:,16], 'k--', label = 'Extracellular Bacteria', linewidth = 3)
plt.semilogy(t, sol[:,19], 'g-.', label = 'Bacterial Leakage', linewidth = 3)
#plt.axis([t[0], t[-1], 10**-2, 10**4])
plt.legend(loc=int(4), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Bacterial Count', fontsize = fontSize)
plt.savefig('bacteria_vs_time_Tcelldep.jpg', dpi=DPI)
plt.show()
    
plt.figure(2)
plt.semilogy(t, sol[:,0], 'c-', label = 'Resting Macrophages', linewidth = 3)
plt.semilogy(t, sol[:,1] - sol[:,2], 'k--', label = 'Infected Macrophages', linewidth = 3)
plt.semilogy(t, sol[:,3], 'g-.', label = 'Activated Macrophages', linewidth = 3)
plt.axis([t[0],t[-1],1,10**6])
plt.legend(loc=int(7), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Macrophage Count', fontsize = fontSize)
plt.savefig('macrophages_vs_time_Tcelldep.jpg', dpi=DPI)
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
plt.savefig('Tcells_vs_time_Tcelldep.jpg', dpi=DPI)
plt.show()

plt.figure(4)
plt.plot(t, sol[:,4], 'c-', label = 'Th0 Cells', linewidth = 3)
#plt.legend(loc=int(4), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('$T_0$ Cell Count', fontsize = fontSize)
plt.savefig('T0_vs_time_Tcelldep.jpg', dpi=DPI)
plt.show()

plt.figure(5)
plt.plot(t, sol[:,10], 'c-', label = r'$F_{\alpha}$', linewidth = 3)
plt.plot(t, sol[:,11], 'k--', label = r'$I_{\gamma}$', linewidth = 3)
plt.plot(t, sol[:,14], 'r--', label = '$I_4$', linewidth = 3)
plt.plot(t, sol[:,13], 'b-', label = '$I_{10}$', linewidth = 3)
plt.plot(t, sol[:,12], 'g-.', label = '$I_{12}$', linewidth = 3)
plt.legend(loc=int(7), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Cytokine Levels', fontsize = fontSize)
plt.axis([t[0],t[-1],0,800])
plt.savefig('cytokines_vs_time_Tcelldep.jpg', dpi=DPI)
plt.show()

plt.figure(6)
plt.semilogy(t, sol[:,17], 'c-', label = 'Collagen', linewidth = 3)
plt.semilogy(t, sol[:,18], 'k--', label = 'MMP', linewidth = 3)
plt.legend(loc=int(4), fontsize = fontSize)
#plt.axis([t[0],t[-1],10,10**4])
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Concentration', fontsize = fontSize)
plt.savefig('Collagen_and_MMP-1_vs_time_Tcelldep.jpg', dpi=DPI)
plt.show()

plt.figure(7)
plt.plot(t, sol[:,19], 'c-', label = 'Bacterial Leakage', linewidth = 3)
#plt.legend(loc=int(7), fontsize = fontSize)
#plt.axis([t[0],t[-1],1,10**4])
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Bacterial Leakage', fontsize = fontSize)
plt.savefig('Bacterial_Leakage_Tcelldep.jpg', dpi=DPI)
plt.show()

plt.figure(8)
plt.plot(t, sol[:,10], 'c-', label = 'TNF', linewidth = 3)
plt.legend(loc=int(4), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Cytokine Levels', fontsize = fontSize)
#plt.axis([t[0],t[-1],0,450])
plt.savefig('TNFalpha_vs_time_Tcelldep.jpg', dpi=DPI)
plt.show()

plt.figure(9)
plt.plot(t, sol[:,17], 'c-', label = 'Collagen', linewidth = 3)
#plt.legend(loc=int(4), fontsize = fontSize)
#plt.axis([t[0],t[-1],10,10**4])
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Collagen Concentration (g$\cdot$cm$^{-3}$)', fontsize = fontSize)
plt.ticklabel_format(style = 'sci', useOffset=True, scilimits = (0,0), axis='y')
#plt.tight_layout()
plt.savefig('Collagen_vs_time_Tcelldep.jpg', dpi=DPI)
plt.show()