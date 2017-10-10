# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 07:47:31 2017

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

def DeleteGene(self,y,t,p):
         
    return 0
    
### Uncomment one of the below ###
### Leaking ###
    
deletedGene = "$M_A$"        
model2.dM_Adt = types.MethodType(DeleteGene, model)

###################################
    
#deletedGene = "$T_8$"        
#model2.dT_8dt = types.MethodType(DeleteGene, model)

###################################
### Latent ###
    
#deletedGene = "$T_2$"        
#model2.dT_2dt = types.MethodType(DeleteGene, model)

###################################
 
#deletedGene = "$I_{10}$"        
#model2.dI_10dt = types.MethodType(DeleteGene, model)

###################################
    
#deletedGene = "$I_4$"        
#model2.dI_4dt = types.MethodType(DeleteGene, model)

###################################
### Active ###
    
#deletedGene = "$T_1$"        
#model2.dT_1dt = types.MethodType(DeleteGene, model)

###################################
    
#deletedGene = "$F_a$"        
#model2.dF_adt = types.MethodType(DeleteGene, model)

###################################
    
#deletedGene = "$I_y$"        
#model2.dI_ydt = types.MethodType(DeleteGene, model)

###################################
### Periodic Switching ##
    
#deletedGene = "$T_0$"        
#model2.dT_0dt = types.MethodType(DeleteGene, model)

###################################
    
#deletedGene = "$T_{80}$"        
#model2.dT_80dt = types.MethodType(DeleteGene, model)

###################################
   
#deletedGene = "$I_{12}$"        
#model2.dI_12dt = types.MethodType(DeleteGene, model)

###################################
###Causes Simulation Errors ###
 
#deletedGene = "$T_c$"        
#model2.dT_cdt = types.MethodType(DeleteGene, model)

###################################


y2,t2,p2,sol2,model2 = RunModel(model=model2, inputMode='file', inputSource='ColDet-model-parameters.csv')


fontSize = 14
DPI = 160

plt.figure(1)
plt.semilogy(t, sol[:,1], 'c-', label = 'Baseline', linewidth = 3)
plt.semilogy(t, sol2[:,1], 'k--', label = (deletedGene + ' Deleted'), linewidth = 3)
#plt.axis([t[0], t[-1], 10**-2, 10**4])
plt.legend(loc=int(4), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Infected Macrophage Count', fontsize = fontSize)
#plt.savefig('M-I-Gene-Deletion.jpg', dpi=DPI)
plt.show()
    
plt.figure(2)
plt.semilogy(t, sol[:,15], 'c-', label = 'Baseline', linewidth = 3)
plt.semilogy(t, sol2[:,15], 'k--', label = (deletedGene + ' Deleted'), linewidth = 3)
#plt.axis([t[0],t[-1],1,10**6])
plt.legend(loc=int(7), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Intracellular Bacteria Count', fontsize = fontSize)
#plt.savefig('B-I-Gene-Deletion.jpg', dpi=DPI)
plt.show()

plt.figure(3)
plt.semilogy(t, sol[:,16], 'c-', label = 'Baseline', linewidth = 3)
plt.semilogy(t, sol2[:,16], 'k--', label = (deletedGene + ' Deleted'), linewidth = 3)
plt.legend(loc=int(4), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Extracellular Bacteria Count', fontsize = fontSize)
#plt.savefig('B-E-Gene-Deletion.jpg', dpi=DPI)
plt.show()

plt.figure(4)
plt.semilogy(t, sol[:,17], 'c-', label = 'Baseline', linewidth = 3)
plt.semilogy(t, sol2[:,17], 'k--', label = (deletedGene + ' Deleted'), linewidth = 3)
plt.legend(loc=int(4), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Collagen', fontsize = fontSize)
#plt.savefig('Collagen-Gene-Deletion.jpg', dpi=DPI)
plt.show()

plt.figure(5)
plt.semilogy(t, sol[:,19], 'c-', label = 'Baseline', linewidth = 3)
plt.semilogy(t, sol2[:,19], 'k--', label = (deletedGene + ' Deleted'), linewidth = 3)
plt.legend(loc=int(4), fontsize = fontSize)
plt.xlabel('Time (Days)', fontsize = fontSize)
plt.ylabel('Bacterial Leakage', fontsize = fontSize)
#plt.savefig('B-L-Gene-Deletion.jpg', dpi=DPI)
plt.show()
