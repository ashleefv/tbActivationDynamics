# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 18:02:24 2017

@author: Steve
"""

from RunModel import RunModel
import matplotlib.pyplot as plt
from WSDetModel import WSDetModel
from DetModel import DetModel
from ColDetModel import ColDetModel
import time

start = time.time()
modelToRun1 = ColDetModel()
y,t,p,solColEqs,model = RunModel(model=modelToRun1, inputMode='file', inputSource='ColDet-model-parameters.csv')
#y,t,p,sol = RunModel(ModelName=WSDetModel, inputMode='file', inputSource='WSDet-model-parameters.csv')
end = time.time()
print(end-start)

start = time.time()
modelToRun2 = ColDetModel(useColEqs=0)
y,t,p,solNoColEqs,model = RunModel(model=modelToRun2, inputMode='file', inputSource='ColDet-model-parameters.csv')
#y,t,p,sol = RunModel(ModelName=WSDetModel, inputMode='file', inputSource='WSDet-model-parameters.csv')
end = time.time()
print(end-start)

fontSize = 18
DPI = 160

plt.figure(1)
plt.semilogy(t, solColEqs[:,15], 'c-', label = 'Intracellular Bacteria', linewidth = 3)
plt.semilogy(t, solColEqs[:,16], 'k--', label = 'Extracellular Bacteria', linewidth = 3)
plt.semilogy(t, solColEqs[:,19], 'g-', label = 'Bacterial Leakage', linewidth = 3)
plt.axis([t[0], t[-1], 10**-2, 10**4])
plt.legend(loc=int(4), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Bacterial Count', fontsize = fontSize)
plt.savefig('bacteria_vs_time_ColEqs.jpg', dpi=DPI)
plt.show()

plt.figure(2)
plt.semilogy(t, solNoColEqs[:,15], 'c-', label = 'Intracellular Bacteria', linewidth = 3)
plt.semilogy(t, solNoColEqs[:,16], 'k--', label = 'Extracellular Bacteria', linewidth = 3)
#plt.semilogy(t, solNoColEqs[:,19], 'g-', label = 'Bacterial Leakage', linewidth = 3)
plt.axis([t[0], t[-1], 10**-2, 10**4])
plt.legend(loc=int(4), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Bacterial Count', fontsize = fontSize)
plt.savefig('bacteria_vs_time_noColEqs.jpg', dpi=DPI)
plt.show()