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


def ParamAdjust(itKey, parameters, y0, t0, t, baselineLeakage, scale, sampleTiming, ModelName=ColDetModel):
    #Run the model with the appropriate changes
    
    parameters[itKey] = parameters[itKey] * scale
    model = ModelName()
    sol = odeint(model.DifEqs, y0, t, args = (parameters,), full_output = 0, atol = 1.0e-1, hmax = t0[2])
    
    while (model.tSwitch > (t[-1] - (sampleTiming*3.5))):
        t = np.arange(t0[0],(t[-1]*1.5)+t0[2],t0[2])
        sol = odeint(model.DifEqs, y0, t, args = (parameters,), full_output = 0, atol = 1.0e-1, hmax = t0[2])
    
    
    switchIndex = math.ceil(model.tSwitch / t0[2])

    indexArray = -1
    
    gain = np.array([((sol[indexArray,19]-baselineLeakage)/baselineLeakage)/(scale-1)])
    #gain is percent change in leakage over percent change in parameter
    
    
    return [itKey, gain]
    
def RunTask(inQueue, outQueue):
    #Function run by the child processes to handle tasks and report results
    while(1):
        func,args = inQueue.get()
        
        outQueue.put(func(*args))
        
        inQueue.task_done()
    
if __name__ == '__main__':
    
    #initialization
    
    numCPUs = multiprocessing.cpu_count()       
    
    inputSource='ColDet-model-parameters.csv'
    scales = (1.05, 1.04, 1.03, 1.02, 1.01, 0.99, 0.98, 0.97, 0.96, 0.95)
    #scales = (1.05,0.95)
    
    #Depreciated
    sampleTiming = 365
    y0, t0, parameters = ImportParamFromFile(inputSource, ColDetModel)

    #Runbase case
    t = np.arange(t0[0],t0[1]+t0[2],t0[2])   
    
    model = ColDetModel()
    sol = odeint(model.DifEqs, y0, t, args = (parameters,), full_output = 0, atol = 1.0e-1, hmax = t0[2])
    
    indexArray = -1    
    baselineLeakage = np.array(sol[indexArray,19])
    
    #reinit parameters to revert changes that may have occured calculation of the model
    y0, t0, parameters = ImportParamFromFile(inputSource, ColDetModel)
    
    #initializing multiprocessing
    processes = []
    resultsArrays= []
    finishQueues = []    
    multiprocessing.freeze_support()    
    runQueue = multiprocessing.JoinableQueue()  
    
    #Set up results Queues and Processes
    for i in range(numCPUs):
        finishQueues.append(multiprocessing.Queue())
        print('Made result Queue #' + repr(i+1))
        processes.append(multiprocessing.Process(target = RunTask, args=(runQueue,finishQueues[i])))
        processes[i].start()
        print('Made process #' + repr(i+1))
        
    #Run sensativity analysis for each scale factor
    paramKeys = list(parameters.keys())        
    
    for j in range(len(scales)):
        
        keyArray = np.empty(0)
        gainArray = np.empty(0)
        
        #Set up and push tasks to queue
        tasks = [(ParamAdjust, (key, parameters, y0, t0, t, baselineLeakage, scales[j], sampleTiming)) for key in paramKeys]
        for task in tasks:
            runQueue.put_nowait(task)
            
        runQueue.join()
        
        #Pull results and append them to a list defined in an outer scope
        results = []
        for i in range(len(finishQueues)):
            while(not finishQueues[i].empty()):
                results.append(finishQueues[i].get())
                
        for i in range(len(results)):
            keyArray = np.append(keyArray, results[i][0])
            gainArray = np.append(gainArray, results[i][1])
        
        #Combine keyArray and gainArray into one array
        compositeArray = np.append(np.transpose(np.atleast_2d(keyArray)), np.transpose(np.atleast_2d(gainArray)), axis=1)
        
        resultsArrays.append(compositeArray)
                

    for i in range(len(processes)):
            processes[i].terminate()  
    #Dump results       
    file = open('LocalSensitivity.p', 'wb')
    pickle.dump(resultsArrays, file)
    file.close()
    
    #Organize results into a spreadsheet        
    resultsOutput = np.empty((len(keyArray)+1, 0))
    dType1 = np.dtype([('param', '<U32'), ('100 Days', '<U32')])
    dType2 = np.dtype('<U32')
    
    for i in range(len(scales)):
        
        resultsArrays[i].dtype = dType1
        resultsArrays[i] = np.sort(resultsArrays[i], axis = 0, order='param')   
        resultsArrays[i].dtype = dType2
        
        dummyArray1 = np.array([['scale', scales[i], '']])
        dummyArray2 = np.full((len(keyArray), 1), '', dtype = '<U5')

        dummyArray2 = np.append(dummyArray2, resultsArrays[i], axis = 1)
        dummyArray3 = np.append(dummyArray1, dummyArray2, axis=0)
        resultsOutput = np.append(resultsOutput, dummyArray3, axis = 1)
    
    np.savetxt('SensativityAnalysis.csv', resultsOutput, fmt = '%s', delimiter = ',')
    