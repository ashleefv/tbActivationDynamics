# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 14:08:58 2017

@author: Steve
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt

file = open('LocalSensitivity.p', 'rb')
resultsArray = pickle.load(file)
file.close()

params = []
for i in range(len(resultsArray[0])):
    params.append(resultsArray[0][i][0])

sensitivities = np.zeros(len(resultsArray[0]))
for i in range(len(resultsArray)):
    for j in range(len(resultsArray[0])):
        sensitivities[j] += float(resultsArray[i][j][1])

sensitivities /= 10

paramsSensitivities = []
for i in range(len(params)):
    paramsSensitivities.append([params[i], sensitivities[i]])
    
paramsSensitivities.sort(key = (lambda entry: np.abs(entry[1])), reverse=True)

params = [x[0] for x in paramsSensitivities]
sensitivities = [x[1] for x in paramsSensitivities]

indexCL = params.index('C_L')
params.pop(indexCL)
sensitivities.pop(indexCL)

indexKC = params.index('k_C') +1

posSensitivities = [np.max([x,0]) for x in sensitivities]

negSensitivities = [np.min([x,0]) for x in sensitivities]

def ParseParam(param):
    retParam = '$'
    bracketFlag = 0
    for char in param:
        if char == 'a':
            retParam += '\\alpha'
            continue
        if char == 'u':
            retParam += '\\mu'
            continue
        if char == '_':
            retParam += '_{'
            bracketFlag = 1
            continue
        retParam += char
    
    if bracketFlag == 1:
        retParam += '}'
    
    retParam += '$'
    
    return retParam

labels = [ParseParam(param) for param in params[0:indexKC]]
            
#params[0] = '$\mu_{TNF}$'
fontSize = 7
DPI = 320
width = 0.5



left = np.arange(1.75,indexKC + 1.75,1)
pos = plt.bar(left, posSensitivities[0:indexKC], width)
neg = plt.bar(left, np.abs(negSensitivities[0:indexKC]), width, color='none')
plt.legend((pos[0], neg[0]), ('Positive Sensititivies', 'Negativie Sensitivities'), fontsize = fontSize)
plt.xticks([x+.25 for x in left], labels, fontsize = fontSize)
plt.yticks(np.arange(0,5.0,0.5), fontsize=fontSize)
plt.ylabel('Absolute Value of Sensitivity', fontsize=fontSize)
plt.xlabel('Parameters', fontsize=fontSize)
plt.savefig('Local_Sensitivities.jpg', dpi=DPI)
plt.show()



