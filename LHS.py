#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 16:41:14 2017

@author: noahgade
"""
import numpy as np
import pandas as pd   

g={"a":.1,"b":2,"c":30,"d":400,"e":5000}

def LHSfunc(nsamp,dictionary):
    
# Set the min and max ranges by factors of the dictionary
    minfactor = 1
    maxfactor = 10
    
    xminlist={}
    xmaxlist={}

    for key, value in dictionary.items():
        xminlist[key] = value * minfactor
        xmaxlist[key] = value * maxfactor

# To overwrite a variable value or range, write here:
    # xminlist["variable"] = value  # for minimum value change
    # xmaxlist["variable"] = value  # for maximum value change

    xminset = pd.Series(xminlist)
    xmin = xminset.values

    xmaxset = pd.Series(xmaxlist)
    xmax = xmaxset.values
    
    nvar = len(xmin)
    
    labels = xminset.index
 
    coeff = (xmax - xmin) / nsamp

    random_array = np.random.random((nsamp,nvar))

    sample_list = np.arange(1,nsamp+1,1)[:, np.newaxis]

    sample_array = np.tile(sample_list,(nvar))

    ordered_array = coeff*(sample_array-random_array) + xmin

    section = np.split(ordered_array,nvar,1)

    for i in range(0,nvar):
        np.random.shuffle(section[i])

    S_list = np.concatenate(section)

    S_finalarray = np.reshape(S_list,(nsamp,nvar),order='F')
    
    varset = S_finalarray.tolist()
              
    output={}
    for i in range(nsamp):
        output[i] = dict(zip(labels,varset[i])) 

    return(output)
    
if __name__ == '__main__':
    print(LHSfunc(9,g))
    




