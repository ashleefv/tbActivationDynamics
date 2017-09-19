# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 11:48:58 2017

@author: Steve
"""



from RunModel import RunModel
import matplotlib.pyplot as plt
from ColDetModel import ColDetModel

y0,t,p,sol1,model1 = RunModel(model=ColDetModel(), inputMode='file', inputSource='ColDet-model-parameters.csv')
#y,t,p,sol = RunModel(ModelName=WSDetModel, inputMode='file', inputSource='WSDet-model-parameters.csv')


p['k_C'] = 2.82E+05
y0,t,p,sol2,model2 = RunModel(model=ColDetModel(), inputMode='arg', paramsIn = [y0, t, p])


p['k_C'] = 5.64E+05
y0,t,p,sol3,model3 = RunModel(model=ColDetModel(), inputMode='arg', paramsIn = [y0, t, p])



fontSize = 14
DPI = 160

plt.figure(1)
plt.plot(t, sol1[:,19], 'c-', label = '$k_C = 1.41 * 10^{5}$', linewidth = 3)
plt.plot(t, sol2[:,19], 'k--', label = '$k_C = 2.82 * 10^{5}$', linewidth = 3)
plt.plot(t, sol3[:,19], 'g-.', label = '$k_C = 5.64 * 10^{5}$', linewidth = 3)
plt.axis([t[0], t[-1],0, 6e3])
plt.legend(loc=int(2), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Bacterial Leakage', fontsize = fontSize)
plt.savefig('leakage_vs_time_kC.jpg', dpi=DPI)
plt.show()

plt.figure(2)
plt.plot(t, sol1[:,16], 'c-', label = '$k_C = 1.41 * 10^{5}$', linewidth = 3)
plt.plot(t, sol2[:,16], 'k--', label = '$k_C = 2.82 * 10^{5}$', linewidth = 3)
plt.plot(t, sol3[:,16], 'g-.', label = '$k_C = 5.64 * 10^{5}$', linewidth = 3)
plt.axis([t[0], t[-1], 0, 1.25e2])
plt.legend(loc=int(4), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Extracellular Bacteria', fontsize = fontSize)
plt.savefig('Be_vs_time_kC.jpg', dpi=DPI)
plt.show()

plt.figure(3)
plt.plot(t, sol1[:,15], 'c-', label = '$k_C = 1.41 * 10^{5}$', linewidth = 3)
plt.plot(t, sol2[:,15], 'k--', label = '$k_C = 2.82 * 10^{5}$', linewidth = 3)
plt.plot(t, sol3[:,15], 'g-.', label = '$k_C = 5.64 * 10^{5}$', linewidth = 3)
plt.axis([t[0], t[-1], 0, 5e3])
plt.legend(loc=int(4), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Intracellular Bacteria', fontsize = fontSize)
plt.savefig('Bi_vs_time_kC.jpg', dpi=DPI)
plt.show()

plt.figure(4)
plt.plot(t, sol1[:,1] - sol1[:,2], 'c-', label = '$k_C = 1.41 * 10^{5}$', linewidth = 3)
plt.plot(t, sol2[:,1] - sol2[:,2], 'k--', label = '$k_C = 2.82 * 10^{5}$', linewidth = 3)
plt.plot(t, sol3[:,1] - sol3[:,2], 'g-.', label = '$k_C = 5.64 * 10^{5}$', linewidth = 3)
plt.axis([t[0], t[-1], 0, 1e2])
plt.legend(loc=int(4), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Infected Macrophages', fontsize = fontSize)
plt.savefig('Mi_vs_time_kC.jpg', dpi=DPI)
plt.show()

plt.figure(5)
plt.plot(t, sol1[:,17], 'c-', label = '$k_C = 1.41 * 10^{5}$', linewidth = 3)
plt.plot(t, sol2[:,17], 'k--', label = '$k_C = 2.82 * 10^{5}$', linewidth = 3)
plt.plot(t, sol3[:,17], 'g-.', label = '$k_C = 5.64 * 10^{5}$', linewidth = 3)
plt.axis([t[0], t[-1], 0, 1e-3])
plt.legend(loc=int(2), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Collagen Concentration', fontsize = fontSize)
plt.savefig('Col_vs_time_kC.jpg', dpi=DPI)
plt.show()

plt.figure(5)
plt.plot(t, sol1[:,18], 'c-', label = '$k_C = 1.41 * 10^{5}$', linewidth = 3)
plt.plot(t, sol2[:,18], 'k--', label = '$k_C = 2.82 * 10^{5}$', linewidth = 3)
plt.plot(t, sol3[:,18], 'g-.', label = '$k_C = 5.64 * 10^{5}$', linewidth = 3)
plt.axis([t[0], t[-1], 0, 2e-8])
plt.legend(loc=int(4), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('MMP Concentration', fontsize = fontSize)
plt.savefig('MMP_vs_time_kC.jpg', dpi=DPI)
plt.show()