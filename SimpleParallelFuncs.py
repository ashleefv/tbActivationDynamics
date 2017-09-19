# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 14:33:56 2017

@author: Steve
"""
#Simple Parallel Pool
#This is an attempt to be a simple wrapper for the multiprocessing pool class
#that uses the parallel map function to run a variety of independent functions
import multiprocessing.pool
import operator

#This is unused.
def _CallFunction(task):
    return task[0](*task[1])
    
#Accepts an iterable where the first element in each element is a func, and the
#second is a tuple with the function arguments.
    
#Ex: [[func1, (arg1,arg2)],
#     [func2, (arg1,)],
#     [func3, (arg1, arg2, arg3)]]
    
def RunTasks(it, numProcessors = None):
    
    multiprocessing.freeze_support()  
    pool =  multiprocessing.pool.Pool(numProcessors)
    
    resultObjects = [pool.apply_async(i[0], i[1]) for i in it]
    
    results = [j.get() for j in resultObjects]
    
    pool.close()
    
    return results



#A test case
if __name__ == "__main__":
    
        
    iteratable = [(operator.add, (1,2)), (operator.sub, (4,3)), (sum, ([5,6,7],))]
    
    results = RunTasks(iteratable, 3)
    
    for i in results:
        print(i)
        
    #input("Press any key to close...")