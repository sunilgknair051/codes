# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 10:43:49 2018

@author: User
"""

#import numpy as np
#import matplotlib.pyplot as plt
#
#n = 20
#X, Y = np.mgrid[0:n, 0:n]
#T = np.arctan2(Y - n / 2., X - n/2.)
#R = 10 + np.sqrt((Y - n / 2.0) ** 2 + (X - n / 2.0) ** 2)
#U, V = R * np.cos(T), R * np.sin(T)
#
#plt.axes([0.025, 0.025, 0.95, 0.95])
#plt.quiver(X, Y, U, V, R, alpha=0.5)
#plt.quiver(X, Y, U, V, edgecolor='k', facecolor='None', linewidth=.5)
#
#plt.xlim(0, n)
#plt.xticks(())
#plt.ylim(0, n)
#plt.yticks(())
#
#plt.show()
#
#



import numpy as np
import matplotlib.pyplot as plt



NSTEPS = 3; NX = 20; NY = NX; ndir = 9; u_max = 0.04; rho0 = 1.; nu = 1./6.; tau = 3.*nu+0.5;
dirx = [0,1,0,-1,0 ,1,-1,-1, 1]
diry = [0,0,1,0 ,-1,1, 1,-1,-1]
w0 = 4./9.
ws = 1./9.
wd = 1./36.
wi = [w0,ws,ws,ws,ws,wd,wd,wd,wd]
f = np.zeros((NX,NY,ndir))
r = np.ones((NX,NY))
u = np.zeros((NX,NY))
v = np.zeros((NX,NY))

n = 20
X, Y = np.mgrid[0:n, 0:n]

#X, Y = np.mgrid[0:NX, 0:NY]
R = 10 + np.sqrt((Y - NY / 2.0) ** 2 + (X - NX / 2.0) ** 2)

plt.axes([0, (2.*np.pi), 0, (2.*np.pi)])
plt.quiver(X, Y, u, v, R, alpha=.5)
plt.quiver(X, Y, u, v, edgecolor='k', facecolor='None', linewidth=.5)

plt.xlim(0, 2.*np.pi)
plt.xticks(())
plt.ylim(0, 2.*np.pi)
plt.yticks(())

plt.show()

      
             