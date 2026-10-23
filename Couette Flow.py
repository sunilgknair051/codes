import numpy as np
import matplotlib.pyplot as plt


def plot(): 
    X, Y = np.meshgrid(np.arange(0,NX), np.arange(0, NY))
    plt.figure(); plt.quiver(Y,X,ux,uy)
    
def init_equilibrium(r,ux,uy):
    f = np.ones((NX,NY,ndir))
    for x in range(NX):
        for y in range(NY):
            for i in range(ndir):
                cidotu = dirx[i]*ux[x,y] + diry[i]*uy[x,y]
                f[x,y,i] = wi[i]*r[x,y]*(1. + 3.*cidotu + 4.5*cidotu**2 - 1.5*(ux[x,y]**2+uy[x,y]**2))
    return f

def stream(f_src):
    f_dst = np.ones((NX,NY,ndir))
    for x in range(NX):
        for y in range(1,NY):
            for i in range(ndir):
                xmd = (NX+x-dirx[i]) % NX
                ymd = (NY+y-diry[i]) % NY
                f_dst[x,y,i] = f_src[xmd,ymd,i]
        y = 0
        for i in [0,1,3,4,7,8]:
            xmd = (NX+x-dirx[i]) % NX
            ymd = (NY+y-diry[i]) % NY
            f_dst[x,y,i] = f_src[xmd,ymd,i]
        for i in [2,5,6]:
            cidotu = dirx[pn[i]]*uLowWallx + diry[pn[i]]*uLowWally 
            f_dst[x,y,i] = f_src[x,y,pn[i]] - 6 * wi[pn[i]]*r[x,y]*cidotu 
#            f_dst[x,y,i] = f_src[x,y,pn[i]] 
        y = NY-1
        for i in [0,1,2,3,5,6]:
            xmd = (NX+x-dirx[i]) % NX
            ymd = (NY+y-diry[i]) % NY
            f_dst[x,y,i] = f_src[xmd,ymd,i]
        for i in [4,7,8]:
            cidotu = dirx[pn[i]]*uUpWallx + diry[pn[i]]*uUpWally 
#            f_dst[x,y,i] = f_src[x,y,pn[i]] - 6 * wi[pn[i]]*r[x,y]*cidotu 
            f_dst[x,y,i] = f_src[x,y,pn[i]]
    return f_dst

def compute_rho_u(f):
    for x in range(NX):
        for y in range(NY):
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
    tauinv = 2./(6.*nu+1.); omtauinv = 1.-tauinv
    for x in range(NX):
        for y in range(NY):
            for i in range(ndir):
                cidotu = dirx[i]*ux[x,y] + diry[i]*uy[x,y]
                feq = wi[i]*r[x,y]*(1. + 3.*cidotu + 4.5*cidotu**2 - 1.5*(ux[x,y]**2+uy[x,y]**2))
                f_out[x,y,i] = omtauinv * f[x,y,i] + tauinv * feq 
                # eq. 3.3
    return f_out
    
scale = 1; NX= 8*scale; NY = NX ; t= 0; ndir = 9; w0 = 4./9.; ws = 1./9.; wd = 1./36.
nu = 1./6.; tau = 3. * nu + 0.5; rho0 = 1;
NSTEPS = 100 * scale**2
wi = [w0,ws,ws,ws,ws,wd,wd,wd,wd]
dirx = [0,1,0,-1, 0,1,-1,-1, 1]
diry = [0,0,1, 0,-1,1, 1,-1,-1]
pn = [0,3,4,1,2,7,8,5,6] 
u_max = 0.04/scale
uLowWallx =  u_max * .5; uLowWally = 0; uUpWallx = 0; uUpWally = 0; 
 
r = np.ones((NX,NY)); ux = 0 * np.ones((NX,NY)); uy = np.zeros((NX,NY))

f1 = init_equilibrium(r,ux,uy)
for n in range(NSTEPS):
    f2 = stream(f1)
    r,ux,uy = compute_rho_u(f2)
    print(f2[0,0,:]-f1[0,0,:],)
    print(np.sum(f2[0,0,:]-f1[0,0,:]))
    print(uy[0,0],r[0,0])
    if n in [0,10,99]:
        plot()
    f1 = collide(f2,r,ux,uy)
    
                    #eq. 5.26, p 180
"""
"""