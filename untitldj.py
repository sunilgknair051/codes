# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 10:02:38 2018

@author: User
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

def FitFunc(time, ExpPars):
    return ExpPars[0] * np.exp(- time / ExpPars[1])

def fit_model(FitParsGuess, x, y):
    def func(FitPars):
        ymodel = FitFunc(x,FitPars)
        return y - ymodel
    FitPars, ier = leastsq(func, FitParsGuess)
    func_fit = FitFunc(time, FitPars)
    return FitPars, func_fit  

time = np.empty(0)
exp_decay = np.empty(0)
for line in open('filll.txt', 'r'):
    values = [float(s) for s in line.split()]
    time = np.append(time, values[0])
    exp_decay = np.append(exp_decay, values[1])
    
ExpParsGuess = np.array([1.,2.]) 
FitPars, exp_decay_fit = fit_model(ExpParsGuess, time, exp_decay) 
plt.plot(time,exp_decay,exp_decay_fit)
print(FitPars)