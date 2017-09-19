# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 12:46:35 2017

@author: Steve
"""

import numpy as np
from scipy.signal import argrelextrema

def periodicConvergenceTest():
    #Develop either some error codes or exceptions
    x = np.array(range(1,601))
    y = np.sin(x)*np.e**-x
    
    
    yMaximaIndicies = np.array(argrelextrema(y, np.greater)).flatten()
    yMaxima = y[yMaximaIndicies]
    xInterval = 12
    tol = 1e-3
    
    yIndiciesDiff = yMaximaIndicies[1:len(yMaximaIndicies)] - yMaximaIndicies[0:-1]
    yMaximaDiff = yMaxima[1:len(yMaxima)] - yMaxima[0:-1]
    
    #Remove sets of maxima grater than xInterval apart
    
    yIndDiffComp = np.greater(yIndiciesDiff, xInterval)
    if np.all(np.equal(yIndDiffComp, np.zeros(len(yIndDiffComp)))):
        #All of the maxima are within xInterval
        lastYIndex = 0
    else:
        lastYIndex = np.nonzero(yIndDiffComp)[0][-1] + 1
        if lastYIndex == len(yMaximaDiff):
            #None of the maxima are within xInterval of each other
            #if errorCode:
                #return (False,1)
            return False
            
            
    yMaximaDiff = yMaximaDiff[lastYIndex:len(yMaximaDiff)]
    
    if np.sum(yMaximaDiff) >= 0:
        return False
    else:
        return True
    
    
if __name__ == '__main__':
    print(periodicConvergenceTest())