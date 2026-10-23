# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 08:45:32 2018

@author: User
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
#import seaborn
#from matplotlib import animation
#import pandas as pd
#import pymc as mc


#LBM in python D2Q9#
#Variables involved#

def taylor_green_i(t,x,y):
      kx = 2.* np.pi/NX
      ky = 2.* np.pi/NY
      td = 1./(nu*(kx**2+ky**2))
      
      X= x+0.5
      Y= y+0.5
      
      ux = -u_max * np.sqrt(ky/kx)*np.cos(kx*X)*np.sin(ky*Y)*np.exp(-1.*t/td)
      uy =  u_max * np.sqrt(kx/ky)*np.sin(kx*X)*np.cos(ky*Y)*np.exp(-1.*t/td)
      
      #Pressure calculation#
      P = -0.25*rho0*u_max*u_max*((ky/kx)*np.cos(2.*kx*X)+(kx/ky)*np.cos(2.*ky*Y))*np.exp(-2.*t/td)
      
      rho = rho0+3.*P
      
      return rho, ux, uy
  
#outer calculation#
def taylor_green_o(t):
          for x in range(NX):
              for y in range(NY):
                 
                   r[x,y],ux[x,y],uy[x,y] =taylor_green_i(t,x,y)
                   
          return r, ux, uy
 
def plot(): 
    X, Y = np.meshgrid(np.arange(0, NX), np.arange(0, NY))
    plt.figure()
    plt.title('Taylor-Green')
    Q = plt.quiver(X, Y, ux, uy)    
    

#Initialization of particle population with equilibrium values#
def init_equilibrium(r, ux, uy):
    for x in range(NX):
        for y in range(NY):
            for i in range(ndir):
                cidotu = dirx[i]*ux[x,y] + diry[i]*uy[x,y]
                f[x,y,i] = wi[i]*r[x,y]*(1. + 3.*cidotu) + 4.5*cidotu**2-1.5*(ux[x,y]**2 + uy[x,y]**2)
    return f            

#Particle Streaming Step#
def stream(f_src):
    for x in range(NX):
        for y in range(NY):
            for i in range(ndir):
                xmd = (NX+x-dirx[i]) % NX
                ymd = (NY+y-diry[i]) % NY
                f_dst[x,y,i] = f_src[xmd,ymd,i];
    return f_dst                

#Computation of density and velocity#
def compute_rho_u(f):
    for x in range(NX):
        for y in range(NY):
            r[x,y] = 0
            ux [x,y] = 0
            uy [x,y] = 0
            for i in range(ndir):
                r[x,y]  = r[x,y] +f[x,y,i]
                ux[x,y] = ux[x,y] + dirx[i] * f[x,y,i]
                uy[x,y] = uy[x,y]  + diry[i] * f[x,y,i]
                ux [x,y] = ux [x,y]/r[x,y]
                uy [x,y] = uy [x,y]/r[x,y]
    return r,ux,uy
                
#Particle Collision Step#
def collide (f, r, ux, uy):
 
    tau_inv = 2./(6.*nu+1.)
    m_tau_inv = 1. - tau_inv
    
    for x in range(NX):
        for y in range(NY):
            for i in range(ndir):
                cidotu   = dirx[i]*ux[x,y] + diry[i]*uy[x,y]
                f_eq     = wi[i]*r[x,y]*(1 + 3*cidotu)+ 4.5*cidotu**2-1.5*(ux[x,y]**2+uy[x,y]**2)
                f_out[x,y,i] = m_tau_inv *f[x,y,i] + tau_inv*f_eq
    return f_out

#Taylor-Green calculation#

#Main function call#    
#def main(t,r,u,v):      
   #for y in range(NY):
      #for x in range(NX):


NSTEPS = 25; scale =1; t = 0; 
NX = 101; NY = NX; ndir = 9;                                     
u_max = 0.02; rho0 = 1.; nu = 1./6.; tau = 3.*nu+0.5;
dirx = [0,1,0,-1,0 ,1,-1,-1, 1]
diry = [0,0,1,0 ,-1,1, 1,-1,-1]
w0 = 4./9.
ws = 1./9.
wd = 1./36.
wi = [w0,ws,ws,ws,ws,wd,wd,wd,wd]
f = np.ones((NX,NY,ndir))
rho = np.ones((NX,NY))
r = np.ones((NX,NY))
ux = np.ones((NX,NY))
uy = np.ones((NX,NY))
f_dst = np.ones((NX,NY,ndir))
f_out = np.ones((NX,NY,ndir))
r, ux, uy = taylor_green_o(0)
plot()
lbm1 = init_equilibrium(r,ux,uy)
f = open('filll.txt', 'wb')#
for n in range(NSTEPS):
    lbm2 = stream(lbm1)
    r,ux,uy = compute_rho_u(lbm2)
    
    print(n,ux[10,15]) #assign to an array
    f.write("%10.4e %10.4e\n" % (n, ux[10,15]))
    lbm1 = collide(lbm2,r,ux,uy)
f.close()      
    
#np.savetxt('filll.txt', zip(n, ux), fmt="%5.2f  %10.4e")  
#np.savetxt('filll.txt', np.column_stack((n,ux[10,15])), fmt='%.18e', delimiter=' ')  
#np.savetxt('filll.txt', fmt='%.18e', delimiter=' ') 
#f = open('filll.txt', 'wb')
#np.savetxt(f, np.column_stack((n,ux[10,15])))
#for i in range(n):
 #   f.write("%10.4e %10.4e\n" % (i, ux[1,2]))
#f.close()    

#file = open("saas.txt", "w")
#for index in range(len(n)):
#    file.write(str(n[index]) + " " + str(ux) + "\n")
#file.close()
#
#
#    
#data = (n, ux[10,15])    
#newdata = np.squeeze(data) # Shape is now: (10, 80)
#plt.plot(data)
    
#np.savetxt('aray3.txt', zip(n, ux[0][0]), fmt="%5.2f  %10.4e")    
#f = open('aray3.txt', 'wb')
#for i in range(len(n)):
#    f.write("%10.4e %10.4e\n" % (n[i], ux[i][i]))
#f.close()    



X, Y = [], []
#for line in open('filll.txt', 'r'):
#  values = [float(s) for s in line.split()]
#  X.append(values[0])
#  Y.append(values[1])
#
#plt.plot(X, Y,label='points_single term',marker='o')
#plt.legend()
#plt.title("Taylor-Green-Plot1")
#plt.xlabel("TSteps")
#plt.ylabel("Values")
##plt.savefig("yes1.pdf")
#plt.show()


#X1, Y1 = [], []
#for line in open('aray3.txt', 'r'):
#  values = [float(s) for s in line.split()]
 # X1.append(values[0])
  #Y1.append(values[1])

#plt.plot(X1, Y1,label='points_added_terms',marker='o')
#plt.legend()
#plt.title("Taylor-Green")
#plt.xlabel("TSteps")
# plt.ylabel("Values")
#plt.savefig("yes2.pdf")
#plt.show()
#
#def plot(x,exp,fit):
#    fig, ax = plt.subplots()
#    ax.plot(x, exp, 'x')
#    ax.plot(x, fit, '-')
#
def FitFunc(time, ExpPars):
    return ExpPars[0] * np.exp(- time / ExpPars[1])

def fit_model(FitParsGuess, x, y): 

    def func(FitPars):
        ymodel = FitFunc(x,FitPars)
        return y - ymodel
    FitPars, ier = leastsq(func, FitParsGuess)
    func_fit = FitFunc(time, FitPars)
    return FitPars, func_fit  

for line in open('filll.txt', 'r'):
  values = [float(s) for s in line.split()]
  time = X.append(values[0])
  exp_decay = Y.append(values[1])
#time      = values[(0)]
#exp_decay = np.array([1,0.8,.7,.6,.55,.47,.43,.39,.35])
ExpParsGuess = np.array([1.,2.]) 
FitPars, exp_decay_fit = fit_model(ExpParsGuess, time, exp_decay) 
plot(time,exp_decay,exp_decay_fit)
print(FitPars)


#
#Define your function
#def func(x, a, b, c):
#    return a*np.exp(-b*x) + c

#def func(time, ExpPars):
#    return ExpPars[0]*np.exp(- time / ExpPars[1])
#
##create a dummy dataset
##x = np.linspace(0,4,50)
#x = np.linspace(0,values[0],10)
##create the 'ideal' negative exponential curve, using the function
##y = func(x, 2.5, 1.3, 0.5)
#y = func(1, 2)
##the dummy data, with some noise thrown in
#yn = y + 0.2*np.random.normal(size=len(x))
#
## Using the python scipy function curve_fit with function and input variables
#popt, pcov = curve_fit(func, x, yn)
#print popt.size
#print pcov.size
#f1=plt.plot(x,yn,x, func(x,popt[0],popt[1],popt[2]))
#plt.title('curve fit')
#plt.xlabel('x-vals')
#plt.ylabel('y-vals')
#plt.show()



