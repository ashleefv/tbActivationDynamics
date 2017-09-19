# Author: Steven Ruggiero
# Email: Stevemruggiero@gmail.com
# 
# The purpose of this file is simply to define a method to import values for parameters
# so further code can be developed assuming that population of a dictionary (or other structure)
# of parameters and their values is handled separately.
#
# I want to rewrite this so it can use an object/structure of some sort to return
# The parameters back in a way specified by the object

import numpy

#Todo: figure out how to put this exception in if not every parameter is defined
class InputError(Exception):
    """Custom exception that is thrown when the inputFile has an error"""
    
    def __init__(self, message):
        self.message = message

def ImportParamFromFile(fileName, model):
# Args:
# fileName: path to the plain text file containing parameters for the model
#
# Returns
# initConditionArray: An array containing all of the initial values for each variable
# timeArray: An array containing the inital and final times and the time step    
# paramDict: A dictionary that accepts a string containing the parameter
#     and returns a float containg the parameter's value
    
    
    fileArray = numpy.loadtxt(fileName, dtype = {'names':('Param', 'Value'), 'formats': ('S10', 'float64')}, delimiter = ',')
    # Opens fileName and loads the contents into a vector of tuples

    
    paramDict = {}
    # Dictionary initialization
    
    # I need to create a more elegant solution than this for arranging the initial
    # conditions into an array
    # Also, this is specifically for the ColDet model currently.
    # Nevermind, just required defining the dict in the class that contains the model
    initConditionDict = model.initConditionDict
    
    initConditionArray = numpy.empty(len(initConditionDict))
    conditionFlags = numpy.zeros(len(initConditionDict))
    
    timeArray = numpy.empty(3)
    timeDict = {'t_0' : 0,
                't_f' : 1,
                't_s' : 2}
    
    inputType = None
    
                  
    def paramInput(lineFromFile):
        param = lineFromFile[0].decode('UTF-8')
        value = lineFromFile[1]
        paramDict[param] = value
        #Populates paramDict with the values from fileArray
        
    def initCondInput(lineFromFile):
        variable = lineFromFile[0].decode('UTF-8')
        index = initConditionDict[variable]
        conditionFlags[index] = 1
        initConditionArray[index] = lineFromFile[1]
        
    def initTimeInput(lineFromFile):
        timeParam = lineFromFile[0].decode('UTF-8')
        index = timeDict[timeParam]
        timeArray[index] = lineFromFile[1]
        
        
    
    inputDict = { 0 : paramInput,
                  1 : initCondInput,
                  2 : initTimeInput}
    # inputDictionary is my way of effectively implementing a switch statement
    
    
    for i in range(0,len(fileArray)):
        if fileArray[i][0].decode('UTF-8') == '_inputType':
            inputType = fileArray[i][1]
            #print(inputType)
        else:
            if inputType == None:
                raise InputError('Input type is not defined')
            else:
                inputFunc = inputDict[inputType]
                inputFunc(fileArray[i])
             
    initCondReturn = numpy.zeros(0)
    for i in range(len(conditionFlags)):
        if conditionFlags[int(i)] == 1:
            initCondReturn = numpy.append(initCondReturn, initConditionArray[int(i)])
    
    return [initCondReturn, timeArray, paramDict]


        
    
