# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 16:53:45 2017

@author: Steve
"""

from ImportParamFromFile import ImportParamFromFile
from LHS import LHSfunc
from ColDetModel import ColDetModel
import math
import numpy as np
from scipy.integrate import odeint
from SimpleParallelFuncs import RunTasks
from SALib.sample import saltelli
from SALib.analyze import sobol
import pickle

def CalcGain(paramDict, y0, t0, t, sampleTiming, ModelName=ColDetModel):
    
    
    model = ModelName()
    sol = odeint(model.DifEqs, y0, t, args = (paramDict,), full_output = 0, atol = 1.0e-1, hmax = t0[2])
    
    #while (model.tSwitch > (t[-1] - (sampleTiming*3.5)) or (model.tSwitch == -1 and not model.colDetFlag)):
    #    t = np.arange(t0[0],(t[-1]*1.5)+t0[2],t0[2])
    #    sol = odeint(model.DifEqs, y0, t, args = (parameters,), full_output = 0, atol = 1.0e-1, hmax = t0[2])
    
    
    #switchIndex = math.ceil(model.tSwitch / t0[2])

    #indexArray = [switchIndex + math.ceil(sampleTiming/t0[2]), switchIndex + math.ceil(2 * sampleTiming/t0[2]), switchIndex + math.ceil(3 * sampleTiming/t0[2])]
    
    bacterialLeakage = sol[-1,18]
    #gain is percent change in leakage over percent change in parameter
    
    return bacterialLeakage
    #return [bacterialLeakage, paramDict]


if __name__ == '__main__':
    numSamples = 25000
    model = ColDetModel()
    inputSource='ColDet-model-parameters.csv'
    #y0, t0, parameters = ImportParamFromFile(inputSource, model)
    
    sampleTiming = 365
    y0, t0, parameters = ImportParamFromFile(inputSource, ColDetModel)
 
    t = np.arange(t0[0],t0[1]+t0[2],t0[2])   
    
    CalcGain(parameters, y0, t0, t, sampleTiming)
#    t = np.array(range(int((t0[1]+t0[2])/t0[2])))
#    t = t0[2] * t
    # t =  range(initialTime, finalTime + timeStep, timeStep)  
    
    #model = ColDetModel()
    #sol = odeint(model.DifEqs, y0, t, args = (parameters,), full_output = 0, atol = 1.0e-1, hmax = t0[2])
    
    #while (model.tSwitch > (t[-1] - (sampleTiming*3.5)) or (model.tSwitch == -1 and not model.colDetFlag)):
    #    t = np.arange(t0[0],(t[-1]*1.5)+t0[2],t0[2])
    #    sol = odeint(model.DifEqs, y0, t, args = (parameters,), full_output = 0, atol = 1.0e-1, hmax = t0[2])
        
    #t = np.arange(t0[0],t0[1]+t0[2],t0[2])
    
    #switchIndex = math.ceil(model.tSwitch / t0[2])
    
    #indexArray = [switchIndex + math.ceil(sampleTiming/t0[2]), switchIndex + math.ceil(2 * sampleTiming/t0[2]), switchIndex + math.ceil(3 * sampleTiming/t0[2])]
    
    #baselineLeakage = np.array(sol[indexArray,18])
    
    globalAnalysisParams = ['u_TNF', 's_4b1', 'u_T1', 'a_30', 'm', 'sr_3bc', 'k_52',
        'sr_3b', 'u_Tc', 'c_52', 'u_iy', 'a_32', 'u_P1',
        'u_i', 's_1', 'k_C']
    
    problem = {'num_vars': len(globalAnalysisParams),
               'names' : globalAnalysisParams,
               'bounds' : [[parameters[globalAnalysisParams[0]]*0.5, parameters[globalAnalysisParams[0]]*1.5 ],
                           [160,200],
                           [parameters[globalAnalysisParams[2]]*0.5, parameters[globalAnalysisParams[2]]*1.5 ],
                           [1e-3,2e-2],
                           [0.5, 1],
                           [1e3, 8e4],
                           [0.07, 1], 
                           [parameters[globalAnalysisParams[7]]*0.5, parameters[globalAnalysisParams[7]]*1.5 ],
                           [parameters[globalAnalysisParams[8]]*0.5, parameters[globalAnalysisParams[8]]*1.5 ],
                           [10, 100],
                           [2.16, 33.2],
                           [parameters[globalAnalysisParams[11]]*0.5, parameters[globalAnalysisParams[11]]*1.5 ],
                           [parameters[globalAnalysisParams[12]]*0.5, parameters[globalAnalysisParams[12]]*1.5 ],
                           [0, 0.005], 
                           [50, 110],
                           [5.206e-2, 1]]}
                           
    param_values = saltelli.sample(problem, numSamples, calc_second_order=True)
    
    paramDicts = []
    for i in range(0, len(param_values)):
        paramDicts.append(dict(parameters))
        for j in range(0,len(param_values[i])):
            paramDicts[i][globalAnalysisParams[j]] = param_values[i][j]
    
    #init Dict to pass to LHSfunc
#    LHSdict = {}
#    for i in globalAnalysisParams:
#        LHSdict[i] = parameters[i]
#    
#    LHSSampleDicts = LHSfunc(numSamples, LHSdict)
#    
#    paramDicts = []
#    
#    for i in range(0, len(LHSSampleDicts)):
#        paramDicts.append(dict(parameters))
#        for key,value in LHSSampleDicts[i].items():
#            paramDicts[i][key] = value
            
    tasks = [(CalcGain, (i, y0, t0, t, sampleTiming)) for i in paramDicts]
    
    results = RunTasks(tasks, numProcessors = 8)
    
    Si = sobol.analyze(problem, np.array(results), print_to_console=False)
        
    file = open('SensitivityResults.p', 'wb')
    pickle.dump(Si, file)
    file.close()

