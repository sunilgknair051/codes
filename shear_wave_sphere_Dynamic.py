import numpy as np; import matplotlib.pyplot as plt; from matplotlib import cm 

def inSample(x,y):  return (x - xSphe)**2 + (y - ySphe)**2 < RadSphe**2
def aboveTop(y):    return y > NY-1
def belowBottom(y): return y < 0
def inLiq(x,y):     return (not inSample(x,y)) and (not aboveTop(y)) and (not belowBottom(y))

def init_equilibrium(dr,ux,uy):
    h = np.zeros((NX,NY,ndir),dtype = complex)
    for x in range(NX):
        for y in range(NY):
            for i in range(ndir):
                cidotu = dirx[i]*ux[x,y] + diry[i]*uy[x,y]
                h[x,y,i] = dr[x,y] + 3.*cidotu
    return h

def stream(h_src):
    h_dst = np.zeros((NX,NY,ndir),dtype = complex)
    FxLiqOnWall = 0; FxLiqOnSphe = 0;
    for x in range(NX):
        for y in range(NY):
            if inLiq(x,y):
                for i in range(ndir):
                    ibar = noslip[i]; xmd = (NX+x-dirx[i]) % NX; ymd = y-diry[i]
                    if (not aboveTop(ymd)) and (not belowBottom(ymd)) and (not inSample(xmd,ymd)) : 
                        h_dst[x,y,i] = h_src[xmd,ymd,i]
                    elif aboveTop(ymd)    : h_dst[x,y,i] = h_src[x,y,ibar]
                    elif belowBottom(ymd) :
                        h_dst[x,y,i] = h_src[x,y,ibar] - 6 * wi[ibar] * dirx[ibar] / wd
                        FxLiqOnWall += (h_dst[x,y,i] + h_src[x,y,ibar])*dirx[i]*wi[i]#*uxWall    # Eqs. 5.79, 5.80
                    elif inSample(xmd,ymd) : 
                        h_dst[x,y,i] = h_src[x,y,ibar] - 6 * wi[ibar] * dirx[ibar] / wd # UWall = 1
                        FxLiqOnSphe += (h_dst[x,y,i] + h_src[x,y,ibar])*dirx[i]*wi[i]
    return h_dst,FxLiqOnWall,FxLiqOnSphe

def compute_drho_u(h_in):
    for x in range(NX):
        for y in range(NY):
            dr[x,y] = 0; ux[x,y] = 0; uy[x,y] = 0;
            for i in range(ndir):
                dr[x,y] += wi[i]*h_in[x,y,i]                   
                ux[x,y] += wi[i]*h_in[x,y,i]*dirx[i]                   
                uy[x,y] += wi[i]*h_in[x,y,i]*diry[i]
    return dr,ux,uy

def collide(h_in,dr,ux,uy):
    h_out = np.ones((NX,NY,ndir),dtype = complex)
    for x in range(NX):
        for y in range(NY):
            for i in range(ndir):
                cidotu = dirx[i]*ux[x,y] + diry[i]*uy[x,y]
                heq = dr[x,y] + 3.*cidotu
                h_out[x,y,i] = omtauinv * h_in[x,y,i] + tauinv * heq + 1j*omega * h_in[x,y,i] - omega2by2 * heq    # eq. 3.3
    return h_out
                    
def Update_uSphe_dfc():
    uxSphe = (SpringConstCont/(1j*omega) + FxLiqOnSphe) / (1j*omega*MassSphe + SpringConstCont/(1j*omega))
    FxSpheOnWall = SpringConstCont/(1j*omega) * (uxSphe - uxWall)
    TanStress = (FxLiqOnWall + FxSpheOnWall) / NX; 
    ZL = TanStress / uxWall; Zq = (rho0/rho_Liq_SI)*(nu/nu_SI)*(delta_SI/delta)*Zq_SI 
    dfc = f0_SI* 1j /np.pi * ZL/Zq
    dfc = -np.conjugate(dfc)
    return dfc,uxSphe,FxSpheOnWall

NX= 20; NY = 20;   
delfac = 1/4.; delta = NY*delfac; k = (1 - 1j)/delta;
nu = 1./6; omega = 2. * nu / delta**2; omega2by2 = 1/2*omega**2; NSTEPS = int(1 / np.abs(omega)); 
RadSphe = NY/6; xSphe = NX/2.; ySphe = 0; uxWall = 1; uxSphe = uxWall; 
MassSphe = 1/2*np.pi*RadSphe**2; SpringConstCont = 10 * omega**2 * MassSphe; 

rho0 = 1.; tau = 3. * nu + 0.5; tauinv = 2./(6.*nu+1.); omtauinv = 1.-tauinv; Zliq = rho0 *(omega*nu)**(1/2);  
nu_SI = 1e-6; delta_SI = 250e-9; rho_Liq_SI = 1e3; Zq_SI = 8.8e6; f0_SI = 5e6; 

ndir = 9; w0 = 4./9.; ws = 1./9.; wd = 1./36.;  
wi = [w0,ws,ws,ws,ws,wd,wd,wd,wd]
dirx = [0,1,0,-1, 0,1,-1,-1, 1]             # 6  2  5      ^ y
diry = [0,0,1, 0,-1,1, 1,-1,-1]             # 3-0-1        | 
noslip =   [0,3,4,1,2,7,8,5,6]              # 7  4  8       -----> x
dr = np.zeros((NX,NY),dtype = complex); 
ux = np.zeros((NX,NY),dtype = complex); 
uy = np.zeros((NX,NY),dtype = complex); 
dfcarr = []

for x in range(NX):
    for y in range(NY): ux[x,y] =  np.exp(-k*y)
h1 = init_equilibrium(dr,ux,uy)

for n in range(NSTEPS):
    h2,FxLiqOnWall,FxLiqOnSphe = stream(h1)
    dr,ux,uy = compute_drho_u(h2)
    h1 = collide(h2,dr,ux,uy)
    dfc,uxSphe,FxSpheOnWall = Update_uSphe_dfc()
    if n%100 == 0 : print(n,"%10.3f"%dfc.real,"%10.3f"%dfc.imag)
    dfcarr = np.append(dfcarr,dfc)

print(NSTEPS,"%10.3f"%dfc.real,"%10.3f"%dfc.imag)
plt.subplot(221); plt.plot(dfcarr[1:].real, 'x'); plt.plot(dfcarr[1:].imag, 'x')
plt.subplot(222); plt.plot(ux[0,:].real, 'x',label='Re(ux)');plt.plot(ux[0,:].imag, 'x',label='Im(ux)')
plt.subplot(223); plt.imshow(ux.real.transpose(),cmap=cm.Reds);
plt.subplot(224); plt.imshow(uy.real.transpose(),cmap=cm.Reds);


