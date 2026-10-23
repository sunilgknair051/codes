# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 09:06:17 2018

@author: User
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
#
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#u, v = np.mgrid[0:2*np.pi:60j, 0:np.pi:160j]
#x = np.cos(u)*np.sin(v)
#y = np.sin(u)*np.sin(v)
#z = np.cos(v)
#
#
#ax.plot_wireframe(x, y, z)

#theta = np.linspace(0, np.pi * 2, 80)
#
#    # equations for great cricles (longitduinal great circles)
#x = R * np.sin(theta[i]) * np.cos((1j / len(theta)) * np.pi * 2)
#y = R * np.sin(theta[i]) * np.sin((1j / len(theta)) * np.pi * 2)
#z = R * np.cos(theta[i])
#
#    # takes the real components of the great circle equation to get the position coords in 3D space
#xr = x.real
#yr = y.real
#


circles = [8., 8., 2.]
#circle=plt.Circle((x, y), r, color='b', fill=False)
#Defining a circle with an array
def get_x_y_co(circles):
    xc = circles[0] #x-co of circle (center)
    yc = circles[1] #y-co of circle (center)
    r = circles[2] #radius of circle
    arr=[]
    for i in range(360):
        y = yc + r*math.cos(i)
        x = xc+ r*math.cos(i)
        x=int(x)
        y=int(y)
        #Create array with all the x-co and y-co of the circle
        arr.append([x,y])
        circles=plt.Circle((x, y), r, color='b', fill=False)
        plt.show()
    return arr


