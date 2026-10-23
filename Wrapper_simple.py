from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import leastsq
import time
import datetime

def Spiral(it, SpiralPars):
    offset = SpiralPars[0] + 1j*SpiralPars[1]
    ampl   = SpiralPars[2] + 1j*SpiralPars[3]
    omc    = SpiralPars[4] + 1j*SpiralPars[5]
    return offset + ampl*np.exp(1j * omc * it)

def fit_model(t_arr, y_arr,SpiralParsGuess):
    y_stacked = np.hstack((np.real(y_arr), np.imag(y_arr)))
    def func(SpiralPars):
        ymodel = Spiral(t_arr,SpiralPars)
        ymodel_stacked = np.hstack((np.real(ymodel), np.imag(ymodel)))
        return y_stacked - ymodel_stacked
    FitPars, ier = leastsq(func,SpiralParsGuess)
    y_arr_fitted = Spiral(t_arr,FitPars)
    chi2 = np.sum((np.abs(y_arr-y_arr_fitted))**2)
    return FitPars,y_arr_fitted,chi2   

def Fit_Spiral():
    f = open('res_ring_in.dat','r')
    lines = f.readlines()
    dfs = np.zeros(len(lines)); dGs = np.zeros(len(lines)); t_arr = np.zeros(len(lines)); 
    for i in range(len(lines)):
        columns = lines[i].split(); t_arr[i] = float(columns[0]); 
        dfs[i] = float(columns[1]); dGs[i] = float(columns[2])
    f.close

    discardini = int(len(dfs)*2./3.)
    t_arr = t_arr[discardini:]-t_arr[discardini] 
    dfs = dfs[discardini:] 
    dGs = dGs[discardini:]
    dfc = dfs + 1j * dGs
    df_fin_Fortr = dfs[-1] 
    dG_fin_Fortr = dGs[-1]

    SpiralParsGuess = np.array([dfs[-1],dGs[-1], dfs[0]-dfs[-1],dGs[0]-dGs[-1],\
                                0.1/len(t_arr),2./len(t_arr)]) 

    SpiralParsFit,dfc_fitted,chi2 = fit_model(t_arr, dfc,SpiralParsGuess)
    df_extrapol = SpiralParsFit[0]
    dG_extrapol = SpiralParsFit[1]
    residual = dfc - dfc_fitted
    if i_plot_Spirals == 1 : 
       plt.subplot(121); plt.plot(t_arr,-dfs,'-'); plt.plot(t_arr,dGs,'-'); 
       plt.plot(t_arr,-dfc_fitted.real,'-'); plt.plot(t_arr,dfc_fitted.imag,'-');
       plt.title(("dfres %7.3f dGamm %7.3f "%(df_extrapol,dG_extrapol)))
       plt.subplot(122); plt.plot(residual.real,'-'); plt.plot(residual.imag,'-'); 
       plt.title(("residuals"))
       plt.show()
    return df_extrapol,dG_extrapol,chi2,df_fin_Fortr,dG_fin_Fortr 

plt.show()

#def Simulate() :
#    f = open('sim_config_B.asc1','w+')
#    f.write("%f\r\n" % target_acc_Hz)
#    f.write("%e\r\n" % rad_Sph_SI)
#    f.write("%e\r\n" % cell_height_SI)
#    f.write("%e\r\n" % grid_reso_SI)
#    f.write("%f\r\n" % cell_diam_by_sphere_diam)
#    f.write("%f\r\n" % truncation_by_diameter)
#    f.write("%d\r\n" % iovt)
#    f.write("%d\r\n" % n_Sph)
#    for i_Sph in range(n_Sph) :
#        f.write("%f %f \r\n" % (xfrac_Sph[i_Sph],zfrac_Sph[i_Sph]))
#    f.close
#    
#    time.sleep(0.1) 
#    f = open('sim_config_B.asc1','w+')
#    f.write("%f\r\n" % target_acc_Hz)
#    f.write("%e\r\n" % rad_Sph_SI)
#    f.write("%e\r\n" % cell_height_SI)
#    f.write("%e\r\n" % grid_reso_SI)
#    f.write("%f\r\n" % cell_diam_by_sphere_diam)
#    f.write("%f\r\n" % truncation_by_diameter)
#    f.write("%d\r\n" % iovt)
#    for i_Sph in range(n_Sph) :
#        f.write("%f %f \r\n" % (xfrac_Sph[i_Sph],zfrac_Sph[i_Sph]))
#    f.close
#    time.sleep(0.1) 
#    
#    os.system("./hello")
#    time.sleep(0.5) 
#    
#    f = open('sim_time_steps.asc1','r')
#    lines = f.readlines()
#    cols = lines[0].split(); nPo = int(cols[2]); nsteps = int(cols[3])
#    cols = lines[1].split(); timesecs = float(cols[0])
#    f.close
#    return nPo,nsteps,timesecs   
#
#def Initialize() :
#    fn_extrapol    = descriptor + datetime.datetime.now().strftime("%m-%d_%H:%M") + '_dfdg_extrapol.asc1'
#    fn_fin_Fortran = descriptor + datetime.datetime.now().strftime("%m-%d_%H:%M") + '_dfdg_fin_Fortran.asc1'
#    fn_header_info = descriptor + datetime.datetime.now().strftime("%m-%d_%H:%M") + '_header_info.asc1'
#
#    f = open(fn_header_info,'w+')
#    f.write("Target Change / Step [Hz] : %f \r\n" %(target_acc_Hz))
#    f.write("Sphere Radius [m] : %e \r\n" %(rad_Sph_SI))
#    f.write("Cell Height [m] : %e \r\n" %(cell_height_SI))
#    f.write("Grid Resolution Ini [m] : %e \r\n" %(grid_reso_SI))
#    f.write("Number of Spheres   : %d \r\n" %(n_Sph))
#    f.write("Cell diameter / Sphere diameter fin : %f \r\n" %(cell_diam_by_sphere_diam))
#    f.write("Truncation / Sphere diameter : %f \r\n" %(truncation_by_diameter))
#    f.close 
#    
#    f = open(fn_extrapol,'w+')
#    f.write('ireso grid_reso ');
#    for i_ovt in range(n_ovt):
#        f.write('OvtOrder df/n[Hz] dG/n[Hz] chi^2 ');
#    f.write("\r\n");     
#    f.close 
#    
#    f = open(fn_fin_Fortran,'w+')
#    f.write('ireso  grid_reso ');
#    for i_ovt in range(n_ovt):
#        f.write('OvtOrder  df/n[Hz] dG/n[Hz] ');
#    f.write("\r\n");     
#    f.close 
#    return fn_extrapol,fn_fin_Fortran
#    
#
#target_acc_Hz = 0.1
#rad_Sph_SI = 6e-9
#cell_height_SI = 10e-9
#grid_reso_SI = 1.2e-9
#truncation_by_diameter = 0.5
#n_Sph = 1
#xfrac_Sph = np.zeros(n_Sph); zfrac_Sph = np.zeros(n_Sph); 
#xfrac_Sph[0] = 0.5
#zfrac_Sph[0] = 0.5
#
i_plot_Spirals = 1
plt.show()
#n_ovt = 1
#
#descriptor = 'simple'
#
#cell_diam_by_sphere_diam = 10
#print'cell_ratio =',cell_diam_by_sphere_diam
#coverage = ((np.pi/4)*(1/(cell_diam_by_sphere_diam**2)))
#print'Coverage =',coverage
#l = [i for i in np.linspace(1,cell_diam_by_sphere_diam,5)]
#
##l = list(np.linspace(1,cell_diam_by_sphere_diam,10))
#print'l=',l
#for a in l :
#
#    fn_extrapol,fn_fin_Fortran = Initialize()
#    
#    f = open(fn_extrapol,'a');    f.write("%i  %e " %(1,cell_diam_by_sphere_diam)); f.close
#    f = open(fn_fin_Fortran,'a'); f.write("%i  %e " %(1,cell_diam_by_sphere_diam)); f.close
#    iovt = 2
#    ovtorder = 2 * iovt -1
#    nPo,nsteps,timesecs = Simulate()
#    df_epol,dG_epol,chi2,df_fin_Fortr,dG_fin_Fortr = Fit_Spiral() 
#    dfbyn = df_epol/ovtorder; dGbyn = dG_epol/ovtorder 
#    dfbyn_fin_Fortr = df_fin_Fortr / ovtorder; dGbyn_fin_Fortr = dG_fin_Fortr / ovtorder 
#    print("%3i, ovt_order %3i steps: %5i Time: %6.3f " %(1,ovtorder,nsteps,timesecs))
#    print("Fitted df/n: %7.3f, dG/n: %7.3f, chisqr : %7.3e" %(dfbyn,dGbyn,chi2))
#    time.sleep(0.1) 
#    f = open(fn_extrapol,'a'); 
#    f.write("%i   %f   %f   %e " %(ovtorder, dfbyn,dGbyn,chi2)); 
#    f.close 
#    time.sleep(0.1) 
#    f = open(fn_fin_Fortran,'a'); 
#    f.write("%i   %f   %f   "    %(ovtorder, dfbyn_fin_Fortr,dGbyn_fin_Fortr)); 
#    f.close 
#    time.sleep(0.1) 
#    f = open(fn_extrapol,'a'); f.write("\r\n"); f.close
#    f = open(fn_fin_Fortran,'a'); f.write("\r\n"); f.close
#

