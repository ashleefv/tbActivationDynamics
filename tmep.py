# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 10:30:00 2017

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

model2=ColDetModel()
model2.baseParameters = dict(p)

#model.parameterAdjustFuncs['a_33'] = lambda y,t,p,pb: 0

def GeneDelete(self,y,t,p):
         
    return 0
        
model2.dT_0dt = types.MethodType(GeneDelete, model)


y2,t2,p2,sol2,model2 = RunModel(model=model2, inputMode='file', inputSource='ColDet-model-parameters.csv')

model3=ColDetModel()
model3.baseParameters = dict(p)

#model.parameterAdjustFuncs['a_33'] = lambda y,t,p,pb: 0

        
model3.dT_8dt = types.MethodType(GeneDelete, model)

y3,t3,p3,sol3,model3 = RunModel(model=model3, inputMode='file', inputSource='ColDet-model-parameters.csv')

model4=ColDetModel()
model4.baseParameters = dict(p)

#model.parameterAdjustFuncs['a_33'] = lambda y,t,p,pb: 0


        
model4.dI_ydt = types.MethodType(GeneDelete, model)

y4,t4,p4,sol4,model4 = RunModel(model=model4, inputMode='file', inputSource='ColDet-model-parameters.csv')


fontSize = 14
DPI = 160

plt.figure(1)
plt.semilogy(t, sol[:,1], 'c-', label = 'Baseline', linewidth = 3)
plt.semilogy(t, sol2[:,1], 'k--', label = '$T_0$ Deleted', linewidth = 3)
plt.semilogy(t, sol3[:,1], 'g-.', label = '$T_8$ Deleted', linewidth = 3)
plt.semilogy(t, sol4[:,1], 'b-', label = '$I_\\gamma$ Deleted', linewidth = 3)
#plt.axis([t[0], t[-1], 10**-2, 10**4])
plt.legend(loc=int(4), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Infected Macrophage Count', fontsize = fontSize)
plt.savefig('M-I-Gene-Deletion.jpg', dpi=DPI)
plt.show()
    
plt.figure(2)
plt.semilogy(t, sol[:,15], 'c-', label = 'Baseline', linewidth = 3)
plt.semilogy(t, sol2[:,15], 'k--', label = '$T_0$ Deleted', linewidth = 3)
plt.semilogy(t, sol3[:,15], 'g-.', label = '$T_8$ Deleted', linewidth = 3)
plt.semilogy(t, sol4[:,15], 'b-', label = '$I_\\gamma$ Deleted', linewidth = 3)
#plt.axis([t[0],t[-1],1,10**6])
plt.legend(loc=int(7), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Intracellular Bacteria Count', fontsize = fontSize)
plt.savefig('B-I-Gene-Deletion.jpg', dpi=DPI)
plt.show()

plt.figure(3)
plt.semilogy(t, sol[:,16], 'c-', label = 'Baseline', linewidth = 3)
plt.semilogy(t, sol2[:,16], 'k--', label = '$T_0$ Deleted', linewidth = 3)
plt.semilogy(t, sol3[:,16], 'g-.', label = '$T_8$ Deleted', linewidth = 3)
plt.semilogy(t, sol4[:,16], 'b-', label = '$I_\\gamma$ Deleted', linewidth = 3)
plt.legend(loc=int(4), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Extracellular Bacteria Count', fontsize = fontSize)
plt.savefig('B-E-Gene-Deletion.jpg', dpi=DPI)
plt.show()

plt.figure(4)
plt.semilogy(t, sol[:,17], 'c-', label = 'Baseline', linewidth = 3)
plt.semilogy(t, sol2[:,17], 'k--', label = '$T_0$ Deleted', linewidth = 3)
plt.semilogy(t, sol3[:,17], 'g-.', label = '$T_8$ Deleted', linewidth = 3)
plt.semilogy(t, sol4[:,17], 'b-', label = '$I_\\gamma$ Deleted', linewidth = 3)
plt.legend(loc=int(4), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Collagen', fontsize = fontSize)
plt.savefig('Collagen-Gene-Deletion.jpg', dpi=DPI)
plt.show()

plt.figure(5)
plt.semilogy(t, sol[:,19], 'c-', label = 'Baseline', linewidth = 3)
plt.semilogy(t, sol2[:,19], 'k--', label = '$T_0$ Deleted', linewidth = 3)
plt.semilogy(t, sol3[:,19], 'g-.', label = '$T_8$ Deleted', linewidth = 3)
plt.semilogy(t, sol4[:,19], 'b-', label = '$I_\\gamma$ Deleted', linewidth = 3)
plt.legend(loc=int(4), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Bacterial Leakage', fontsize = fontSize)
plt.savefig('B-L-Gene-Deletion.jpg', dpi=DPI)
plt.show()
