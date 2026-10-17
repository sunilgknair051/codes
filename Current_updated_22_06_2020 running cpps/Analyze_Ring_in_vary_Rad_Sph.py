import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import os
import configparser

def Get_Labels() : 
    keys =  ["ncycles_max", "target_acc_Hz", "rad_Sph_SI", "cellheight_SI", \
	"gridreso_SI", 	"cell_diam_by_sphere_diam", "frac_above_surf", \
	"nufac", "tandelinv", "iovt", "nx", "nz", "ncycles_max", "outInterval","t_computation_secs"];
    return keys

def Read_Header_Info(filename,keys):
    config = configparser.ConfigParser(); config.read(filename); sectionname = header_filename;
    Header_Info = dict(zip(keys,np.zeros(len(keys))))
    for key in keys:                                                                              
        if key in ['ncycles_max','iovt', 'nx', 'nz', 'ncycles_max', 'outInterval']:
            Header_Info[key] = int(config.get(sectionname, key,fallback=-1))
        else :
            Header_Info[key] = float(config.get(sectionname, key,fallback=-9999.))
    return Header_Info
	
def Write_fit_Results():
    fname = os.getcwd() + '\\' + 'results/' + header_filename + '_dfcs.txt'
    f = open(fname,'w+'); f.close();
    f = open(fname,'a+'); 
    for i in range(len(pars_1)):
        for j in range(len(pars_2)):
            print('{:d} {:d} {:.2e} {:.2e} {:.2f}  {:.2f}'.\
              format(i, j, pars_1[i], pars_2[j], \
                     dfc_arr_12[i,j].real, dfc_arr_12[i,j].imag));
            f.write('{:d} {:d} {:.2e} {:.2e} {:.2f} \t {:.2f} \n'.\
              format(i, j, pars_1[i], pars_2[j], \
                     dfc_arr_12[i,j].real, dfc_arr_12[i,j].imag))      
    f.close(); 

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

def Read_Displacement(filename):
    if os.path.isfile(filename) :
        data = np.loadtxt(filename)
        nx = np.int(np.max(data[:,0]) + 1)
        ny = np.int(np.max(data[:,1]) + 1)
    else :    
        nx = 2
        ny = 2
    ux = np.zeros((nx,ny),dtype = complex);
    uy = np.zeros((nx,ny),dtype = complex);
    uz = np.zeros((nx,ny),dtype = complex);
    if os.path.isfile(filename) :
        for i in range(len(data)):
            x = np.int(data[i,0])
            y = np.int(data[i,1])
            ux[x,y] = data[i,2] + 1j * data[i,3] 
            uy[x,y] = data[i,4] + 1j * data[i,5] 
            uz[x,y] = data[i,6] + 1j * data[i,7] 
    return ux,uy,uz;    

def Fit_Spiral(filename):
    data = np.loadtxt(filename); dfs = data[:,0]; dGs = data[:,1]; t_arr = np.arange(len(dfs));
    discardini = int(len(dfs)*1./3.)
    t_arr = t_arr[discardini:] - t_arr[discardini] 
    dfs = dfs[discardini:] 
    dGs = dGs[discardini:]
    dfc = dfs + 1j * dGs
    SpiralParsGuess = np.array([dfs[-1],dGs[-1], dfs[0]-dfs[-1],dGs[0]-dGs[-1],\
                                0.1/len(t_arr),2./len(t_arr)]) 

    SpiralParsFit,dfc_fitted,chi2 = fit_model(t_arr, dfc,SpiralParsGuess)
    df_extrapol = SpiralParsFit[0]
    dG_extrapol = SpiralParsFit[1]
    residual = dfc - dfc_fitted
    if Show_Plots and (n_sims % gap_display == 0): 
        plt.rcParams["figure.figsize"] = (7,3.5)
        plt.subplot(221); plt.plot(t_arr,-dfs,'-'); plt.plot(t_arr,dGs,'-'); 
        plt.plot(t_arr,-dfc_fitted.real,'-'); plt.plot(t_arr,dfc_fitted.imag,'-');
        plt.title(("df/n:%7.1f  dG/n:%7.1f "%(df_extrapol,dG_extrapol)))
        plt.subplot(222); plt.plot(residual.real,'-'); plt.plot(residual.imag,'-'); 
        plt.title(("residuals"))
        plt.tight_layout();
        plt.show()
    return df_extrapol,dG_extrapol,chi2

keys = Get_Labels()
Show_Plots = True;
header_filename = open('header_filename.txt', 'r').read();
gap_display = 1;

dfc_arr = []; comp_time_arr = []
chi2_arr = []; Header_Info_arr = []; 
n_sims = 0; i_file_exists = 1;
while ((i_file_exists == 1) and n_sims < 1000):
    HInfo_filename    = os.getcwd() + '\\' + 'results/' + header_filename + '_HeaderInfo_' + str(n_sims).zfill(3) + '.dat';
    dfc_by_n_filename = os.getcwd() + '\\' + 'results/' + header_filename + '_ring_in_'    + str(n_sims).zfill(3) + '.dat';
    Displacement_filename     = os.getcwd() + '\\' + 'results/' + header_filename + '_Displ_'       + str(n_sims).zfill(3) + '.dat';

    if os.path.isfile(HInfo_filename) :
        ux,uy,uz                      = Read_Displacement(Displacement_filename);
        x1 = np.arange(len(ux[0,:]))
        y1 = np.arange(len(ux[:,0]))
        X, Y = np.meshgrid(x1, y1)
#        plt.subplot(221); cp = plt.contourf(Y,X,np.transpose(ux.real)); plt.colorbar(cp); 
#        plt.subplot(222); cp = plt.contourf(Y,X,np.transpose(ux.imag)); plt.colorbar(cp); 
#        plt.subplot(223); cp = plt.contourf(Y,X,np.transpose(uy.real)); plt.colorbar(cp); 
#        plt.subplot(224); cp = plt.contourf(Y,X,np.transpose(uy.imag)); plt.colorbar(cp); 
        plt.subplot(221); cp = plt.contourf(Y,X,ux.real); plt.colorbar(cp); 
        plt.subplot(222); cp = plt.contourf(Y,X,ux.imag); plt.colorbar(cp); 
        plt.subplot(223); cp = plt.contourf(Y,X,uy.real); plt.colorbar(cp); 
        plt.subplot(224); cp = plt.contourf(Y,X,uy.imag); plt.colorbar(cp); 
        plt.tight_layout();
        plt.show();

        Header_Info                   = Read_Header_Info(HInfo_filename,keys)
        comp_time = np.float(Header_Info['t_computation_secs'])
        df_extrapol,dG_extrapol,chi2  = Fit_Spiral(dfc_by_n_filename);
        Header_Info_arr = np.append(Header_Info_arr, [Header_Info],                    axis=0);    
        dfc_arr         = np.append(dfc_arr,         [df_extrapol + 1j * dG_extrapol], axis=0);  
        chi2_arr        = np.append(chi2_arr,        [chi2],                           axis=0); 
        comp_time_arr   = np.append(comp_time_arr,   [comp_time],                      axis=0); 
        n_sims += 1;
        print('frac_above_surf',Header_Info['frac_above_surf'])
        print('nufac',Header_Info['nufac'])
    else : i_file_exists = 0;
if n_sims > 1 : 
    for key in keys:
        if Header_Info_arr[1][key] != Header_Info_arr[0][key] and \
          key != 't_computation_secs': key_1 = key

n_dimensions = 1    
for key in keys:
    pars = np.zeros(n_sims)
    for i_sim in range(n_sims) : pars[i_sim] = Header_Info_arr[i_sim][key]
    if (np.std(pars) > 1e-12) and key != key_1 and key != 't_computation_secs' and key != 'nx': 
        n_dimensions = 2; 
        key_2 = key;

print('n_dimensions',n_dimensions)        

if n_dimensions == 1: 
    pars_1 = np.zeros(n_sims)
    for i_sim in range(n_sims) : pars_1[i_sim] = Header_Info_arr[i_sim][key_1]
    print('i ' + key_1 + ' df dG')
    f = open(header_filename + '_dfcs.txt','w+'); f.close();
    f = open(header_filename + '_dfcs.txt','a+'); 
    for i in range(len(pars_1)):
        print('{:d} {:.2e} {:.2f} {:.2f} '.\
              format(i, pars_1[i], dfc_arr[i].real, dfc_arr[i].imag));
        f.write('{:d} {:.2e} \t {:.2f} \t {:.2f}  \n' .\
              format(i, pars_1[i], dfc_arr[i].real, dfc_arr[i].imag));      
    plt.subplot(131); plt.plot(pars_1,dfc_arr.real,'-'); 
    plt.xlabel(key_1)
    plt.title('df/n ' + key_1)
    plt.subplot(132); plt.plot(pars_1,dfc_arr.imag,'-'); 
    plt.xlabel(key_1)
    plt.title('dG/n {:.1e}  {:.1e}'.format(pars_1[0],pars_1[-1]) )
    plt.subplot(133); plt.plot(pars_1,comp_time_arr,'-'); 
    plt.title('comp_time [s]')
    plt.xlabel(key_1)

    plt.show();
    
if n_dimensions == 2: 
    y_start = Header_Info_arr[0][key_2]
    i_sim = 0;
    while i_sim < n_sims-1 : 
        i_sim += 1;
        diff = np.abs(Header_Info_arr[i_sim][key_1] - Header_Info_arr[0][key_1])
        if diff < 1e-12 : break
    n1 = i_sim;
    n2 = np.int(n_sims/n1);
    pars_1 = np.zeros(n1); pars_2 = np.zeros(n2)
    for i1 in range(n1) : pars_1[i1] = Header_Info_arr[i1]    [key_1]
    for i2 in range(n2) : pars_2[i2] = Header_Info_arr[i2*n1] [key_2]
    dfc_arr_12 = np.zeros((n1,n2),dtype = complex)
    for i_sim in range(n1 * n2):
        i1 = i_sim % n1;
        i2 = np.int(i_sim/n1);
        dfc_arr_12[i1,i2] = dfc_arr[i_sim]

   
    plt.subplot(121);
    for j in range(len(pars_1)):
        plt.plot(pars_2,dfc_arr_12[     j,:].real,'-',color = 'b'); 
    plt.xlabel(key_2)
    plt.title('df/n ' + key_1)

    plt.subplot(122);
    for j in range(len(pars_1)):
        plt.plot(pars_2,dfc_arr_12[     j,:].imag,'-',color = 'b'); 
    plt.xlabel(key_2)
    
    gridreso_SI = Header_Info['gridreso_SI']
    plt.title('dG/n {:.1e}  {:.1e} gr {:.1f}'.format(pars_1[0],pars_1[-1],gridreso_SI*1e9 ))
    plt.tight_layout();
    plt.show()
         
    plt.subplot(121);
    for j in range(len(pars_2)):
        plt.plot(pars_1[:-1],dfc_arr_12[:-1,     j].real,'-',color = 'b'); 
    plt.xlabel(key_1)
    plt.title('df/n ' + key_2)

    plt.subplot(122);
    for j in range(len(pars_2)):
        plt.plot(pars_1[:-1],dfc_arr_12[:-1,     j].imag,'-',color = 'b'); 
    plt.xlabel(key_1)
    plt.title('dG/n {:.1e}  {:.1e}'.format(pars_2[0],pars_2[-1]) )
    plt.tight_layout();
    plt.show()
#    Write_fit_Results()
#    dfc_arr_12_liq_sub = np.zeros((len(dfc_arr_12)+1,len(dfc_arr_12[0])+1),dtype = complex)
#    for i in range(len(dfc_arr_12)):
#        dfc_arr_12_liq_sub[i+1,0] = pars_1[i]
#        for j in range(len(dfc_arr_12[0])):
#            dfc_arr_12_liq_sub[i+1,j+1] = dfc_arr_12[i,j] - 716.8 * (-1 +1j)
#            dfc_arr_12_liq_sub[0,j+1] = pars_2[j] 
#    print(dfc_arr_12_liq_sub)
    
         