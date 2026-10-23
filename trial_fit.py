# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 10:45:20 2018

@author: User
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def func(x, a, b, c):
    return ExpPars[0] * np.exp(- time / ExpPars[1])

#Define the data to be fit with some noise:

xdata = np.linspace(0, 10, 100)
y = func(xdata, 2.5, 1.3, 0.5)
np.random.seed(1729)
y_noise = 0.2 * np.random.normal(size=xdata.size)
ydata = y + y_noise
plt.plot(xdata, ydata, 'b-', label='data')

#Fit for the parameters a, b, c of the function func:

popt, pcov = curve_fit(func, xdata, ydata)
popt 
#array([ 2.55423706,  1.35190947,  0.47450618])
plt.plot(xdata, func(xdata, *popt), 'r-',
          label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))

#Constrain the optimization to the region of 0 <= a <= 3, 0 <= b <= 1 and 0 <= c <= 0.5:

popt, pcov = curve_fit(func, xdata, ydata, bounds=(0, [3., 1., 0.5]))
popt 
#array([ 2.43708906,  1.        ,  0.35015434])
plt.plot(xdata, func(xdata, *popt), 'g--',
         label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))



plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()

