import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq




def plot(x,exp,fit):
    fig, ax = plt.subplots()
    ax.plot(x, exp, 'x')
    ax.plot(x, fit, '-')

def FitFunc(time, ExpPars):
    return ExpPars[0] * np.exp(- time / ExpPars[1])

def fit_model(FitParsGuess, x, y): 

    def func(FitPars):
        ymodel = FitFunc(x,FitPars)
        return y - ymodel
    FitPars, ier = leastsq(func, FitParsGuess)
    func_fit = FitFunc(time, FitPars)
    return FitPars, func_fit  
 
time      = np.array([1,2  ,3  ,4,5 ,6,  7,  8,  9])
exp_decay = np.array([1,0.8,.7,.6,.55,.47,.43,.39,.35])
ExpParsGuess = np.array([1.,2.]) 
FitPars, exp_decay_fit = fit_model(ExpParsGuess, time, exp_decay) 
plot(time,exp_decay,exp_decay_fit)
print(FitPars)