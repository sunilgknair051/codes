from __future__ import division
import numpy as np
#import matplotlib.pyplot as plt

import matplotlib.pyplot as plt; from matplotlib import cm


def arrows_plot(): 
    X, Y = np.meshgrid(np.arange(0, NX), np.arange(0,NY))
    plt.figure(1)
    plt.title('Sphere in Liquid')
    plt.quiver(Y, X, ux, uy)
#   change to DPx, DPy plt.quiver(Y, X, ux, uy)
    plt.quiver(Y, X, DPx, DPy)
    
def init_equilibrium(r,ux,uy):
    f = np.zeros((NX,NY,ndir))
    for x in range(NX):
        for y in range(NY):
            if inliquid(x,y):
                for i in range(ndir):
                    cidotu = dirx[i]*ux[x,y] + diry[i]*uy[x,y]
                    if i_linearize == 0:
                        f[x,y,i] = wi[i]*r[x,y]*(1. + 3.*cidotu + 4.5*cidotu**2 - 1.5*(ux[x,y]**2+uy[x,y]**2))
                        #print f
                    elif i_linearize == 1:     
                        f[x,y,i] = wi[i]*r[x,y]*(1. + 3.*cidotu) # run till x^1 to check decay
    return f

def inliquid(x,y): return ((x-x_sphere)**2 + (y-y_sphere)**2 > rad_sphere**2 )
                         
def stream(f_src):
    f_dst = np.zeros((NX,NY,ndir))
    for x in range(NX):
        for y in range(NY):
            for i in range(ndir):
                xmd = (NX+x-dirx[i]) % NX; ymd = (NY+y-diry[i]) % NY
                if inliquid(x,y):
                    if inliquid(xmd, ymd): 
                        f_dst[x,y,i] = f_src[xmd,ymd,i]
                    else :
                        ibar = pn[i]
                        f_dst[x,y,i] = f_src[x,y,ibar] - 6 * wi[ibar]*(dirx[ibar]*ux_sphere + diry[ibar]*uy_sphere)
                        # calculate the f(x,y,i)  following eq 5.79  
                        ci = dirx[i] + diry[i]
                        for x in range(NX):
                            for y in range(NY):
                                for i in range(ndir):
                                    f(x,y,i) += (DPx+DPy)*ci/Dx**3
                        # Dx Dt are one for now
    return f_dst

def forces():
    forces (x,y,i)

def compute_rho_u(f):
    for x in range(NX):
        for y in range(NY):
            if inliquid(x,y):
                r[x,y] = 0; ux[x,y] = 0; uy[x,y] = 0;
                for i in range(ndir):
                    r[x,y]  += f[x,y,i]                   
                    ux[x,y] += f[x,y,i]*dirx[i]                   
                    uy[x,y] += f[x,y,i]*diry[i]
                ux[x,y] = ux[x,y]/r[x,y]                    
                uy[x,y] = uy[x,y]/r[x,y]
    return r,ux,uy

def collide(f,r,ux,uy):
    f_out = np.ones((NX,NY,ndir))
    tauinv = 2./(6.*nu+1.)
    omtauinv = 1.-tauinv
    for x in range(NX):
        for y in range(NY):
            if inliquid(x,y):
               for i in range(ndir):
                  cidotu = dirx[i]*ux[x,y] + diry[i]*uy[x,y]
                  if i_linearize == 0:
                     feq = wi[i]*r[x,y]*(1. + 3.*cidotu + 4.5*cidotu**2 - 1.5*(ux[x,y]**2+uy[x,y]**2))
                  elif i_linearize == 1:     
                    feq = wi[i]*r[x,y]*(1. + 3.*cidotu)
                    f_out[x,y,i] = omtauinv * f[x,y,i] + tauinv * feq
    return f_out
                    
scale = 1; NX= 16*scale; NY = NX; ndir = 9; w0 = 4./9.; ws = 1./9.; wd = 1./36.; i_linearize = 1; NSTEPS = 100* scale**2
wi = [w0,ws,ws,ws,ws,wd,wd,wd,wd]
dirx = [0,1,0,-1, 0,1,-1,-1, 1]         # 6  2  5      ^ y
diry = [0,0,1, 0,-1,1, 1,-1,-1]         # 3-0-1        | 
pn =   [0,3,4,1,2,7,8,5,6]              # 7  4  8       -----> x
nu = 1./6.; tau = 3. * nu + 0.5; rho0 = 1.; u_max = 0.01/scale




#Sphere centre, velocity and radius values #preferred
x_sphere = NX/2; y_sphere = NY/2; ux_sphere = 10.; uy_sphere = 0. ; rad_sphere = 2.
Dx = 1.; Dt = 1.; #Equal spacings
r = np.ones((NX,NY)); ux = np.zeros((NX,NY)); uy = np.zeros((NX,NY))
# define momentum transfers DPx ad DPy like velocities here as zeros
DPx = np.zeros((NX,NY)); DPy = np.zeros((NX,NY));

f1 = init_equilibrium(r,ux,uy)
for n in range(NSTEPS):
    f2 = stream(f1)
    r,ux,uy = compute_rho_u(f2)
    f1 = collide(f2,r,ux,uy)
# calculate total force as sum over all x,y,i
for x in range(NX):
    for y in range(NY):
        for i in range(ndir):
                FxLBM += f(x,y,i)
    
arrows_plot()


eta = nu*rho0; FxorceAnalytical = 4.*np.pi*eta*rad_sphere*ux_sphere;
print FxorceAnalytical
print FxLBM

dFro = 20.*np.pi
print (dFro)
print (dFro-FxorceAnalytical)/FxorceAnalytical


##Plotting a circle using given circle dimensions
circle=plt.Circle((x_sphere, y_sphere), rad_sphere, color='b', fill=False)
plt.axis('scaled')
fig = plt.gcf()
fig.gca().add_artist(circle)
plt.show()
 
plt.clf(); plt.imshow(np.sqrt(ux**2+uy**2).transpose(),cmap=cm.Reds)





#def momentum():   #Calculation of momenta   
#     for x in range(NX):
#        for y in range(NY):
#    
#    
#    
#def particle_forces():  #Particle force calculation
#


