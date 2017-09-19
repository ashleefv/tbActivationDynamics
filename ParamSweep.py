# Author: Steven Ruggiero
# Email: Stevemruggiero@gmail.com

from ImportParamFromFile import ImportParamFromFile
from scipy.integrate import odeint
from ColDetModel import ColDetModel
from RunModel import RunModel
import multiprocessing
import math
import numpy as np
import sys


def ParamSweep(ModelName=ColDetModel, inputMode='file', inputSource='ColDet-model-parameters.csv', scale = 1.05):
    if (inputMode == 'file'):
        y0, t0, parameters = ImportParamFromFile(inputSource)
 
    t = np.arange(t0[0],t0[1]+t0[2],t0[2])   
#    t = np.array(range(int((t0[1]+t0[2])/t0[2])))
#    t = t0[2] * t
    # t =  range(initialTime, finalTime + timeStep, timeStep)  
    
    
    sol = odeint(ModelName, y0, t, args = (parameters,), full_output = 0, atol = 1.0e-1, hmax = t0[2])
    
    while (parameters['tSwitch'] > (t[-1] - 350)):
        t = np.arange(t0[0],(t[-1]*1.5)+t0[2],t0[2])
        sol = odeint(ModelName, y0, t, args = (parameters,), full_output = 0, atol = 1.0e-1, hmax = t0[2])
        
    t = np.arange(t0[0],t0[1]+t0[2],t0[2])
        
    switchIndex = math.ceil(parameters['tSwitch'] / t0[2])
    
    indexArray = [switchIndex + math.ceil(100/t0[2]), switchIndex + math.ceil(200/t0[2]), switchIndex + math.ceil(300/t0[2])]
    
    baselineLeakage = np.array(sol[indexArray,18])
    
    #reinit parameters to revert changes that may have occured calculation of the model
    if (inputMode == 'file'):
        y0, t0, parameters = ImportParamFromFile(inputSource)
    
    keyArray = np.empty(0)
    gainArray = np.empty([0,3])
    
    print(scale)    
    
    for itKey in parameters:
        if (itKey == 'colFlag'):
            continue
        loopParameters = dict(parameters)
        loopParameters[itKey] = loopParameters[itKey] * scale
        
        sol = odeint(ModelName, y0, t, args = (loopParameters,), full_output = 0, atol = 1.0e-1, hmax = t0[2])
        
        while (loopParameters['tSwitch'] > (t[-1] - 350)):
            t = np.arange(t0[0],(t[-1]*1.5)+t0[2],t0[2])
            sol = odeint(ModelName, y0, t, args = (loopParameters,), full_output = 0, atol = 1.0e-1, hmax = t0[2])
        
        t = np.arange(t0[0],t0[1]+t0[2],t0[2])
        
        switchIndex = math.ceil(loopParameters['tSwitch'] / t0[2])
    
        indexArray = [switchIndex + math.ceil(100/t0[2]), switchIndex + math.ceil(200/t0[2]), switchIndex + math.ceil(300/t0[2])]
        
        keyArray =  np.append(keyArray, itKey)
        gain = np.array([((sol[indexArray,18]-baselineLeakage)/baselineLeakage)/((1-scale)/scale)])
        #gain is percent change in leakage over percent change in parameter
        gainArray = np.append(gainArray, gain, axis=0)
    
    
    return [keyArray, gainArray]
    
    
    
