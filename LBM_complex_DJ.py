import numpy as np; import matplotlib.pyplot as plt; from matplotlib import cm 

def inSample(x,y)  : return (x-x_sphere)**2 + (y-y_sphere)**2 < rad_sphere**2
def aboveTop(y)    : return y > NY - 1
def belowBottom(y) : return y < 0
def inLiquid(x,y)  : return (not aboveTop(y)) and (not belowBottom(y)) and (not inSample(x,y))

def init_equilibrium(dr,ux,uy):
    h = np.zeros((NX,NY,ndir),dtype = complex)
    for x in range(NX):
        for y in range(NY):
            for i in range(ndir):
                cidotu = dirx[i]*ux[x,y] + diry[i]*uy[x,y]
                h[x,y,i] = wi[i] * dr[x,y]*(1. + 3.*cidotu) 
    return h
                         
def stream(h_src):
    h_dst = np.zeros((NX,NY,ndir),dtype = complex)
    FxLiqOnWall = 0; FxLiqOnSphe = 0;
    for x in range(NX):
        for y in range(NY):
            for i in range(ndir):
                if inLiquid(x,y):
                    ibar = noslip[i]; xmd = (NX+x-dirx[i]) % NX; ymd = y-diry[i]
                    if inLiquid(xmd, ymd): h_dst[x,y,i] = h_src[xmd,ymd,i]
                    elif aboveTop(ymd)   : h_dst[x,y,i] = h_src[x,y,ibar]
                    elif belowBottom(ymd):
                        h_dst[x,y,i] = h_src[x,y,ibar] - 6 * wi[ibar]*dirx[ibar]
                        FxLiqOnWall += (h_dst[x,y,i] + h_src[x,y,ibar])*dirx[i]
                    elif inSample(xmd,ymd) : 
                        h_dst[x,y,i] = h_src[x,y,ibar] - 6 * wi[ibar] * dirx[ibar]
                        FxLiqOnSphe += (h_dst[x,y,i] + h_src[x,y,ibar])*dirx[i]
    return h_dst,FxLiqOnWall,FxLiqOnSphe

def compute_dr_u(h_in):
    for x in range(NX):
        for y in range(NY):
            if inLiquid(x,y):
                dr[x,y] = 0; ux[x,y] = 0; uy[x,y] = 0;
                for i in range(ndir):
                    dr[x,y] += h_in[x,y,i]                   
                    ux[x,y] += h_in[x,y,i]*dirx[i]                   
                    uy[x,y] += h_in[x,y,i]*diry[i]
    return dr,ux,uy

def collide(h_in,dr,ux,uy):
    h_out = np.ones((NX,NY,ndir),dtype = complex)
    for x in range(NX):
        for y in range(NY):
            if inLiquid(x,y):
                for i in range(ndir):
                    cidotu = dirx[i]*ux[x,y] + diry[i]*uy[x,y]
                    heq = wi[i]*(dr[x,y] + 3.*cidotu)   
                    h_out[x,y,i] = omtauinv * h_in[x,y,i] + tauinv * heq + 1j*omega * h_in[x,y,i] - omega**2/2 * heq # eq. 3.3
    return h_out

def calc_dfc():
    TanStress = (FxLiqOnWall + FxLiqOnSphe)/ NX;
    ZL = TanStress / uxWall; 
    Zq = (rho0/rho_Liq_SI)*(nu/nu_SI)*(delta_SI/delta)*Zq_SI 
    dfc = f0_SI* 1j /np.pi * ZL/Zq
    dfc = -np.conjugate(dfc)
    return dfc
            
NY = 20; ndir = 9; w0 = 4./9.; ws = 1./9.; wd = 1./36.; 
y_sphere = 0; rad_sphere = NY/5.
rad_sphere = NY/10.
 
wi = [w0,ws,ws,ws,ws,wd,wd,wd,wd]
dirx = [0,1,0,-1, 0,1,-1,-1, 1]         # 6  2  5        ^ y
diry = [0,0,1, 0,-1,1, 1,-1,-1]         # 3- 0 -1        | 
noslip =   [0,3,4,1,2,7,8,5,6]          # 7  4  8        |-----> x
nu = 1./6.; tau = 3. * nu + 0.5; 
 

delfac = 1/6.; delta = NY*delfac; omega = 2.*nu / delta**2; 
tauinv = 2./(6.*nu+1.); omtauinv = 1.-tauinv; uxWall = 1;

Zq_SI = 8.8e6; nu_SI = 1e-6; rho0 = 1.; rho_Liq_SI = 1000.; f0_SI = 5.e6; delta_SI = 250.e-9

NSTEPS = int(3./omega)
dfcarr = []
for NX in [20,30,40,50]:
    x_sphere = NX/2;   
    dr = np.zeros((NX,NY),dtype = complex); ux = np.zeros((NX,NY),dtype = complex); uy = np.zeros((NX,NY),dtype = complex);
    h1 = init_equilibrium(dr,ux,uy)
    for n in range(NSTEPS):
        h2,FxLiqOnWall,FxLiqOnSphe = stream(h1)
        dr,ux,uy = compute_dr_u(h2)
        h1 = collide(h2,dr,ux,uy)
        dfc = calc_dfc()
        dfcarr = np.append(dfcarr,dfc)

    print(NSTEPS,"%10.3f"%dfc.real,"%10.3f"%dfc.imag)
    plt.subplot(221); plt.plot(dfcarr[1:].real, 'x'); plt.plot(dfcarr[1:].imag, 'x')
    plt.subplot(222); plt.plot(ux[0,:].real, 'x',label='Re(ux)');plt.plot(ux[0,:].imag, 'x',label='Im(ux)')
    plt.subplot(223); plt.imshow(ux.real.transpose(),cmap=cm.Reds);
    plt.subplot(224); plt.imshow(uy.real.transpose(),cmap=cm.Reds);




