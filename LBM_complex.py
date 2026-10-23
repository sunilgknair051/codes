import numpy as np
import matplotlib.pyplot as plt; from matplotlib import cm


def init_equilibrium(r,ux,uy):
    h = np.zeros((NX,NY,ndir))
    for x in range(NX):
        for y in range(NY):
            if inliquid(x,y):
                for i in range(ndir):
                    cidotu = dirx[i]*ux[x,y] + diry[i]*uy[x,y]
                    if i_linearize == 0:
                        h[x,y,i] = wi[i]*dr[x,y]*(1. + 3.*cidotu + 4.5*cidotu**2 - 1.5*(ux[x,y]**2+uy[x,y]**2))
                    elif i_linearize == 1:     
                        h[x,y,i] = wi[i]*dr[x,y]*(1. + 3.*cidotu) # run till x^1 to check decay
    return h

def inliquid(x,y): return ((x-x_sphere)**2 + (y-y_sphere)**2 > rad_sphere**2 )
                         
def stream(h_src):
    h_dst = np.zeros((NX,NY,ndir),dtype = complex)
    for x in range(NX):
        for y in range(NY):
            DPx[x,y] = 0; DPy[x,y] = 0;  
            for i in range(ndir):
                xmd = (NX+x-dirx[i]) % NX; ymd = (NY+y-diry[i]) % NY
                if inliquid(x,y):
                    if inliquid(xmd, ymd): 
                        h_dst[x,y,i] = h_src[xmd,ymd,i]
                    else :
                        ibar = pn[i]
                        h_dst[x,y,i] = h_src[x,y,ibar] - 6 * wi[ibar]*(dirx[ibar]*ux_sphere + diry[ibar]*uy_sphere)
                        DPx[x,y] += (h_src[x,y,ibar] + h_dst[x,y,i]) * dirx[i]
                        DPy[x,y] += (h_src[x,y,ibar] + h_dst[x,y,i]) * diry[i]
    return h_dst,DPx,DPy

def compute_rho_u(h_in):
    for x in range(NX):
        for y in range(NY):
            if inliquid(x,y):
                dr[x,y] = 0; ux[x,y] = 0; uy[x,y] = 0;
                for i in range(ndir):
                    dr[x,y]  += h_in[x,y,i]                   
                    ux[x,y] += h_in[x,y,i]*dirx[i]                   
                    uy[x,y] += h_in[x,y,i]*diry[i]
                ux[x,y] = ux[x,y]/dr[x,y]                    
                uy[x,y] = uy[x,y]/dr[x,y]
    return dr,ux,uy

def collide(h_in,dr,ux,uy):
    h_out = np.ones((NX,NY,ndir))
    tauinv = 2./(6.*nu+1.)
    omtauinv = 1.-tauinv
    for x in range(NX):
        for y in range(NY):
            if inliquid(x,y):
               for i in range(ndir):
                  cidotu = dirx[i]*ux[x,y] + diry[i]*uy[x,y]
                  if i_linearize == 0:
                     heq = dr[x,y] + 3.*cidotu
                  elif i_linearize == 1:     
                    heq = wi[i]*dr[x,y]*(1. + 3.*cidotu)
                    h_out[x,y,i] = omtauinv * h_in[x,y,i] + tauinv * heq + 1j*omega * h_in[x,y,i] - omega**2 * heq # eq. 3.3
    #+uy_of_pointy[n] = uy[point_x,point_y]
    return h_out
            

            
            
scale = 1; NX= 64*scale; NY = NX; ndir = 9; w0 = 4./9.; ws = 1./9.; wd = 1./36.; i_linearize = 1; 
NSTEPS = 10* scale**2
       
#ux_of_pointx = np.ones(NSTEPS);uy_of_pointy = np.ones(NSTEPS)
#point_x = 16; point_y = 5;   
 
wi = [w0,ws,ws,ws,ws,wd,wd,wd,wd]
dirx = [0,1,0,-1, 0,1,-1,-1, 1]         # 6  2  5      ^ y
diry = [0,0,1, 0,-1,1, 1,-1,-1]         # 3-0-1        | 
pn =   [0,3,4,1,2,7,8,5,6]              # 7  4  8       -----> x
nu = 1./6.; tau = 3. * nu + 0.5; rho0 = 1.; 



x_sphere = NX/2.; y_sphere = NY/2.; ux_sphere = 5.; uy_sphere = 0. ; rad_sphere = NX/3. 
dr = np.zeros((NX,NY),dtype = complex); ux = np.zeros((NX,NY),dtype = complex); uy = np.zeros((NX,NY),dtype = complex);
DPx = np.zeros((NX,NY)); DPy = np.zeros((NX,NY));

delfac = 1/3.; delta = NY*delfac; 
omega = 2.*nu / delta**2; 

h1 = init_equilibrium(dr,ux,uy)
Fx_of_n = np.zeros(NSTEPS); Fy_of_n = np.zeros(NSTEPS)  
ux_of_n = np.zeros(NSTEPS); 
for n in range(NSTEPS):
    h2,DPx,DPy = stream(h1)
    dr,ux,uy = compute_rho_u(h2)
    h1 = collide(h2,dr,ux,uy)
    FxLBM = 0; FyLBM = 0; uxsum = 0; norm = 0 # calculate total force as sum over all x,y,i
    for x in range(NX):
        for y in range(NY):
            FxLBM += DPx[x,y]; FyLBM += DPy[x,y]; 
            if inliquid(x,y) : uxsum += ux[x,y]; norm += 1;
            
    ux_of_n[n] = uxsum / norm; Fx_of_n[n] = FxLBM; Fy_of_n[n] = FyLBM
   #print (np.sqrt(((point_x-x_sphere)**2+(point_y-y_sphere)**2))-rad_sphere)
    
    #x_sphere += ux_sphere
    #ux_of_pointx += ux_sphere
    #print (x_sphere, y_sphere)
plt.clf(); plt.imshow(ux.transpose(),cmap=cm.Reds)
fig, ax = plt.subplots(); ax.plot(Fx_of_n[2:-1], 'k-')

#print ux_of_pointx;
#print uy_of_pointy
