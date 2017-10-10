# Author: Steven Ruggiero
# Email: Stevemruggiero@gmail.com
# 
# This file contains the function that handles collecting input,
# calling, and handleing output from scipy.integrate.odeint.
# Odeint uses lsoda from the FORTRAN odepack library.
#
#
# Arguments:
# ModelName: A function handle that points to the Model to evaluate
# inputMode: A string specifying how the input is being fed to the model
# inputSource: The source of the input. Currently can only be a string with a
#   filename

from ImportParamFromFile import ImportParamFromFile
from scipy.integrate import odeint
from ColDetModel import ColDetModel
import numpy as np
def RunModel(model=ColDetModel(), inputMode='file', inputSource='Deterministic-model-parameters.csv', paramsIn=None):

    
    
    if (inputMode == 'file'):
        y0, t0, parameters = ImportParamFromFile(inputSource, model)
        t = np.arange(t0[0],t0[1]+t0[2],t0[2])  
        
    if (inputMode == 'arg'):
        y0, t, parameters = paramsIn
     
    sol = odeint(model.DifEqs, y0, t, args = (parameters,), full_output = 0, atol = 1.0e-1)
    

    
    return [y0, t, parameters, sol, model]