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
import pickle
import time


def ParamAdjust(itKey, parameters, y0, t0, t, baselineLeakage, scale, ModelName=ColDetModel):
    
    
     
    
    
        
    parameters[itKey] = parameters[itKey] * scale
    
    sol = odeint(ModelName, y0, t, args = (parameters,), full_output = 0, atol = 1.0e-1, hmax = t0[2])
    
    while (parameters['tSwitch'] > (t[-1] - 350)):
        t = np.arange(t0[0],(t[-1]*1.5)+t0[2],t0[2])
        sol = odeint(ModelName, y0, t, args = (parameters,), full_output = 0, atol = 1.0e-1, hmax = t0[2])
    
    
    switchIndex = math.ceil(parameters['tSwitch'] / t0[2])

    indexArray = [switchIndex + math.ceil(100/t0[2]), switchIndex + math.ceil(200/t0[2]), switchIndex + math.ceil(300/t0[2])]
    
    gain = np.array([((sol[indexArray,18]-baselineLeakage)/baselineLeakage)/((1-scale)/scale)])
    #gain is percent change in leakage over percent change in parameter
    
    
    return [itKey, gain]
    
def RunTask(inQueue, outQueue):
    while(not inQueue.empty()):
        func,args = inQueue.get()
        #print(func)
        #print(args)
        outQueue.put(func(*args))
        #For some reason, it seems some weirdness happens with pushing to finishQueue. This sleep command is here to help mitigate that
        time.sleep(0.1)
        inQueue.task_done()
    
if __name__ == '__main__':
    
    inputSource='ColDet-model-parameters.csv'
    scale = 1.05
    
    y0, t0, parameters = ImportParamFromFile(inputSource)
 
    t = np.arange(t0[0],t0[1]+t0[2],t0[2])   
#    t = np.array(range(int((t0[1]+t0[2])/t0[2])))
#    t = t0[2] * t
    # t =  range(initialTime, finalTime + timeStep, timeStep)  
    
    
    sol = odeint(ColDetModel, y0, t, args = (parameters,), full_output = 0, atol = 1.0e-1, hmax = t0[2])
    
    while (parameters['tSwitch'] > (t[-1] - 350)):
        t = np.arange(t0[0],(t[-1]*1.5)+t0[2],t0[2])
        sol = odeint(ColDetModel, y0, t, args = (parameters,), full_output = 0, atol = 1.0e-1, hmax = t0[2])
        
    t = np.arange(t0[0],t0[1]+t0[2],t0[2])
        
    switchIndex = math.ceil(parameters['tSwitch'] / t0[2])
    
    indexArray = [switchIndex + math.ceil(100/t0[2]), switchIndex + math.ceil(200/t0[2]), switchIndex + math.ceil(300/t0[2])]
    
    baselineLeakage = np.array(sol[indexArray,18])
    
    #reinit parameters to revert changes that may have occured calculation of the model
    y0, t0, parameters = ImportParamFromFile(inputSource)
    
    keyArray = np.empty(0)
    gainArray = np.empty([0,3])
    
    multiprocessing.freeze_support()
    
    numCPUs = multiprocessing.cpu_count()
    #Forcing numCPUs to be less than max so I can use my computer while this runs.
    #numCPUs = 4
    runQueue = multiprocessing.JoinableQueue()
    
    
    paramKeys = list(parameters.keys())
    
    paramKeys.remove('colFlag')
    
    tasks = [(ParamAdjust, (key, parameters, y0, t0, t, baselineLeakage, scale)) for key in paramKeys]
    for task in tasks:
        runQueue.put_nowait(task)
        
    finishQueues = []
    
    processes = []
    for i in range(numCPUs):
        finishQueues.append(multiprocessing.Queue())
        print('Made result Queue #' + repr(i+1))
        processes.append(multiprocessing.Process(target = RunTask, args=(runQueue,finishQueues[i])))
        processes[i].start()
        print('Made process #' + repr(i+1))
    
    runQueue.join()
    
    for i in range(len(processes)):
        processes[i].terminate()    
    
    results = []
    for i in range(len(finishQueues)):
        while(not finishQueues[i].empty()):
            results.append(finishQueues[i].get())
            
    for i in range(len(results)):
        keyArray = np.append(keyArray, results[i][0])
        gainArray = np.append(gainArray, results[i][1], axis=0)
        
    resultArray = np.append(np.transpose(np.atleast_2d(keyArray)), gainArray, axis=1)
    
        
