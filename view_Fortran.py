import matplotlib.pyplot as plt; 
from matplotlib import cm 
import numpy as np

np.
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




