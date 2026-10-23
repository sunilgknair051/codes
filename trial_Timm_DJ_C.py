import numpy as np
import matplotlib.pyplot as plt
#from scipy.optimize import leastsq

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


NSTEPS = 100; scale =1; t = 0; 
NX = 32; NY = NX; ndir = 9;                                     
u_max = 0.001 ; rho0 = 1.; nu = 1./6.; tau = 3.*nu+0.5;
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
f = open('fill1.txt', 'wb')#
for n in range(NSTEPS):
    lbm2 = stream(lbm1)
    r,ux,uy = compute_rho_u(lbm2)
    
    print(n,ux[10,15]) #assign to an array
    f.write("%10.4e %10.4e\n" % (n, ux[10,15]))
    lbm1 = collide(lbm2,r,ux,uy)
f.close()      
    


#X, Y = [], []

