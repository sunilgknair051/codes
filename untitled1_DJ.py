# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 10:15:27 2018

@author: User
"""
import numpy as np
import matplotlib.pyplot as plt

NSTEPS = 3; NX = 32; NY = NX; ndir = 9; u_max = 0.04; rho0 = 1.; nu = 1./6.; tau = 3.*nu+0.5;
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

def taylor_green_inner(t,x,y):
      kx = (2.* np.pi)/NX
      ky = (2.* np.pi)/NY
      td = 1./(nu*(kx**2+ky**2))
      
      X= x+0.5
      Y= y+0.5
      
      ux = -u_max * np.sqrt(ky/kx)*np.cos(kx*X)*np.sin(ky*Y)*np.exp(-1.*t/td)
      uy =  u_max * np.sqrt(kx/ky)*np.sin(kx*X)*np.cos(ky*Y)*np.exp(-1.*t/td)
      
      #Pressure calculation#
      P = -0.25*rho0*u_max*u_max*((ky/kx)*np.cos(2.*kx*X)+(kx/ky)*np.cos(2.*ky*Y))*np.exp(-2.*t/td)
      
      rho = rho0+3.*P
      
      r = rho
      u = ux
      v = uy 
      return r,u,v
      
      
def taylor_green_outer(t):
    for y in range(NY):
         for x in range(NX):
             r,u,v = taylor_green_inner(t,x,y)
    return r,u,v   

r,ux,uy = taylor_green_outer(0)

X, Y = np.meshgrid(np.arange(0, NX), np.arange(0, NY))
plt.figure()
plt.title('hello')
Q = plt.quiver(X, Y, ux, uy)
"""
X, Y = np.mgrid[0:NX, 0:NY]
R = 10 + np.sqrt((Y - NY / 2.0) ** 2 + (X - NX / 2.0) ** 2)

plt.axes([0, 2*np.pi, 0, 2*np.pi])
plt.quiver(X, Y, u, v, R, alpha=.5)
plt.quiver(X, Y, u, v, edgecolor='k', facecolor='None', linewidth=.5)

plt.xlim(-1, NX)
plt.xticks(())
plt.ylim(-1, NY)
plt.yticks(())

plt.show()

"""      
             