import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import os
import configparser
Show_Plots = True;
Colors = ['k', 'r', 'g', 'b', 'm', 'y', 'k', 'r','g','b','m','y'];
def Get_Labels() : 
    keys =  ["ncycles_max", "target_acc_Hz", "rad_Sph_SI", "cellheight_SI", \
	"gridreso_SI", 	"cell_diam_by_sphere_diam", "frac_above_surf", \
	"nufac", "tandelinv", "iovt", "nx", "nz", "ncycles_max", "outInterval","jqmax","t_computation_secs"];
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

def Read_Displacement_Sph(filename):
    data = np.transpose(np.loadtxt(filename))
    gridreso_SI = np.float(Header_Info['gridreso_SI']);
    rho_SI = 1e3;
    Mass_per_site = rho_SI * gridreso_SI**3;  #The concept here that we don't have the whole sphere and use the cube instead
    if len(data) > 0 :
        x = data[0];
        y = data[1];
        z = data[2];
        ux = data[3] + 1j * data[4]
        uy = data[5] + 1j * data[6]
        uz = data[7] + 1j * data[8]
        
        x_avg_CM = np.average(x);
        y_avg_CM = np.average(y); 
        z_avg_CM = np.average(z);
        
       #print(x_avg_CM,y_avg_CM,z_avg_CM)
        Mom_Inertia_SI = Mass_per_site * np.sum((x-x_avg_CM)**2 + (y-y_avg_CM)**2) * gridreso_SI**2;
       #print('Mom_Inertia_SI: ',Mom_Inertia_SI)
        # Translation
        ux_avg = np.average(ux);
        uy_avg = np.average(uy);
        uz_avg = np.average(uz);
        # Angular momentum
        L_x = Mass_per_site * (uy * (z - z_avg_CM) - uz * (y - y_avg_CM)) * gridreso_SI**2
        L_y = Mass_per_site * (uz * (x - x_avg_CM) - ux * (z - z_avg_CM)) * gridreso_SI**2
        L_z = Mass_per_site * (ux * (y - y_avg_CM) - uy * (x - x_avg_CM)) * gridreso_SI**2
        # Rotation
        Om_x = np.sum(L_x) / Mom_Inertia_SI 
        Om_y = np.sum(L_y) / Mom_Inertia_SI 
        Om_z = np.sum(L_z) / Mom_Inertia_SI 
    else :    
        ux_avg = 1 + 0*1j;
        uy_avg = 0 + 0*1j;
        uz_avg = 0 + 0*1j;
        Om_x = 0 + 0*1j;
        Om_y = 0 + 0*1j;
        Om_z = 0 + 0*1j;
    return ux_avg, uy_avg, uz_avg, Om_x, Om_y, Om_z     



keys = Get_Labels()
header_filename = 'test013';
n_sims_max = 186;
HInfo_filename             = os.getcwd() + '\\' + 'results/' + header_filename +  '_HeaderInfo_' + str(n_sims_max).zfill(3) + '.dat';
Displacement_Sph_filename  = os.getcwd() + '\\' + 'results/' + header_filename + '_Displ_Sph_'   + str(n_sims_max).zfill(3) + '.dat';
Header_Info                   = Read_Header_Info(HInfo_filename,keys)
ux_avg, uy_avg, uz_avg, Om_x, Om_y, Om_z = Read_Displacement_Sph(Displacement_Sph_filename);
     
#print(' ux = {:.3f} + i {:.3f} \n uy = {:.3f} + i {:.3f} \n uz = {:.3f} + i {:.3f}' \
#      .format(ux_avg.real,ux_avg.imag,uy_avg.real,uy_avg.imag,uz_avg.real,uz_avg.imag))
#print(' Om_x = {:.2e} + i {:.2e} \n Om_y = {:.2e} + i {:.2e} \n Om_z = {:.2e} + i {:.2e}' \
#      .format(Om_x.real,Om_x.imag,Om_y.real,Om_y.imag,Om_z.real,Om_z.imag))
print(Header_Info)
	
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
    data = np.loadtxt(filename); 
    print(data[0])
    dfs = data[:,0]; dGs = data[:,1]; t_arr = np.arange(len(dfs));
    discardini = int(len(dfs)*1./3.)
    t_arr = t_arr[discardini:] - t_arr[discardini] 
    dfs = dfs[discardini:] 
    dGs = dGs[discardini:]
    dfc = dfs + 1j * dGs
    if len(dfs) > 8 :
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
    else : df_extrapol,dG_extrapol,chi2 = np.nan,np.nan,np.nan       
    return df_extrapol,dG_extrapol,chi2

keys = Get_Labels()

#header_filename = open('header_filename.txt', 'r').read();
gap_display = 1;
 
dfc_arr = []; comp_time_arr = []
chi2_arr = []; Header_Info_arr = [];
ux_avg_arr = []; Om_z_arr = [];
n_sims = 0; i_file_exists = 1;
#while ((i_file_exists == 1) and n_sims < 1000):
while ((i_file_exists == 1) and n_sims < n_sims_max):
    HInfo_filename    = os.getcwd() + '\\' + 'results/' + header_filename + '_HeaderInfo_' + str(n_sims).zfill(3) + '.dat';
    dfc_by_n_filename     = os.getcwd() + '\\' + 'results/' + header_filename + '_ring_in_'    + str(n_sims).zfill(3) + '.dat';
    dfc_by_n_ref_filename = os.getcwd() + '\\' + 'results/' + header_filename + '_ring_in_'    + str(n_sims).zfill(3) + '_ref.dat';
    Displacement_filename     = os.getcwd() + '\\' + 'results/' + header_filename + '_Displ_'       + str(n_sims).zfill(3) + '.dat';
    Displacement_Sph_filename = os.getcwd() + '\\' + 'results/' + header_filename + '_Displ_Sph_'   + str(n_sims).zfill(3) + '.dat';

    if os.path.isfile(HInfo_filename) :
        if Show_Plots : 
            ux,uy,uz                      = Read_Displacement(Displacement_filename);
            x1 = np.arange(len(ux[0,:]))
            y1 = np.arange(len(ux[:,0]))
            X, Y = np.meshgrid(x1, y1)
            plt.subplot(221); cp = plt.contourf(Y,X,ux.real); plt.colorbar(cp); 
            plt.subplot(222); cp = plt.contourf(Y,X,ux.imag); plt.colorbar(cp); 
            plt.subplot(223); cp = plt.contourf(Y,X,uy.real); plt.colorbar(cp); 
            plt.subplot(224); cp = plt.contourf(Y,X,uy.imag); plt.colorbar(cp); 
            plt.tight_layout();
            plt.show();

        Header_Info                   = Read_Header_Info(HInfo_filename,keys)
        comp_time = np.float(Header_Info['t_computation_secs'])
        df_extrapol,dG_extrapol,chi2  = Fit_Spiral(dfc_by_n_filename);
        if os.path.isfile(dfc_by_n_ref_filename) : 
            df_extrapol_ref,dG_extrapol_ref,chi2  = Fit_Spiral(dfc_by_n_ref_filename);
        else : df_extrapol_ref,dG_extrapol_ref = df_extrapol,dG_extrapol    
        print("cell_diam_by_sphere_diam",Header_Info["cell_diam_by_sphere_diam"]);
        print("iovt",Header_Info["iovt"]);
#        iovt = Header_Info["iovt"];
#        ovtorder = 2 * iovt + 1; 
#        df_extrapol_ref = -716.8 / ovtorder**0.5
#        dG_extrapol_ref =  716.8 / ovtorder**0.5
        print('df_extrapol,df_extrapol_ref',df_extrapol,df_extrapol_ref)
        print('dG_extrapol,dG_extrapol_ref',dG_extrapol,dG_extrapol_ref)
        df_extrapol -= df_extrapol_ref;
        dG_extrapol -= dG_extrapol_ref;
        if os.path.isfile(Displacement_Sph_filename) :
            ux_avg, uy_avg, uz_avg, Om_x, Om_y, Om_z = Read_Displacement_Sph(Displacement_Sph_filename)
        
        
        Header_Info_arr = np.append(Header_Info_arr, [Header_Info],                    axis=0);    
        dfc_arr         = np.append(dfc_arr,         [df_extrapol + 1j * dG_extrapol], axis=0);  
        chi2_arr        = np.append(chi2_arr,        [chi2],                           axis=0); 
        comp_time_arr   = np.append(comp_time_arr,   [comp_time],                      axis=0); 
        ux_avg_arr      = np.append(ux_avg_arr,      [ux_avg],                         axis=0); 
        Om_z_arr        = np.append(Om_z_arr,        [Om_z],                           axis=0); 
        n_sims += 1;
    #   print('frac_above_surf',Header_Info['frac_above_surf'])
     #  print('nufac',Header_Info['nufac'])
    else : i_file_exists = 0;

if n_sims > 1 : 
    for key in keys:
        if Header_Info_arr[1][key] != Header_Info_arr[0][key] and \
          key != 't_computation_secs' and key != 'nx' and key != 'nz' : 
              key_1 = key
              print('key_1',key_1,Header_Info_arr[1][key],Header_Info_arr[0][key])
n_dimensions = 1    

for key in keys:
    pars = np.zeros(n_sims)
    for i_sim in range(n_sims) : pars[i_sim] = Header_Info_arr[i_sim][key]
    if (np.std(pars) > 1e-12) and key != key_1 and key != 't_computation_secs' and key != 'nx' and key != 'nz' and key != 'jqmax': 
        n_dimensions = 2; 
        key_2 = key;
        print('pars',pars)

print('n_dimensions',n_dimensions)        

if n_dimensions == 1: 

    pars_1 = np.zeros(n_sims)
    for i_sim in range(n_sims) : pars_1[i_sim] = Header_Info_arr[i_sim][key_1]
    print('i ' + key_1 + ' df dG')
    f = open(header_filename + '_dfcs.txt','w+'); f.close();
    f = open(header_filename + '_dfcs.txt','a+'); 
    for i in range(len(pars_1)):
        print('{:d} {:.5e} {:.5f} {:.5f} '.\
              format(i, pars_1[i], dfc_arr[i].real, dfc_arr[i].imag));
        f.write('{:d} {:.5e} \t {:.5f} \t {:.5f}  \n' .\
              format(i, pars_1[i], dfc_arr[i].real, dfc_arr[i].imag));      

    #when n_dimensions = 1
    if key_1 == "cell_diam_by_sphere_diam":
        for i_sim in range(n_sims) : 
            rad_Sph_SI = Header_Info["rad_Sph_SI"];
            gridreso_SI = Header_Info["gridreso_SI"];
            nx = Header_Info_arr[i_sim]["nx"];
            nz = Header_Info_arr[i_sim]["nz"];
            key_1 = "coverage"
            coverage = np.pi* (rad_Sph_SI / gridreso_SI)**2 / (nx * nz)
            pars_1[i_sim] = coverage
        key_1 = "coverage"
#         coverage = pars_1**-2 * np.pi / 4.; 
#         pars_1 = coverage

    plt.subplot(321);
    for j in range(len(pars_1)): plt.plot(pars_1,dfc_arr.real,'-',color = 'b'); 
    plt.title('$\Delta$f/n ' + key_1)
    
    plt.subplot(322);
    for j in range(len(pars_1)): plt.plot(pars_1,dfc_arr.imag,'-',color = 'b'); 
    plt.tight_layout();
    
    plt.subplot(323);
    for j in range(len(pars_1)): plt.plot(pars_1,ux_avg_arr.real - 1,'-',color = 'b'); 
    plt.title('ux_avg - 1')
    plt.tight_layout();

    plt.subplot(324);
    for j in range(len(pars_1)): plt.plot(pars_1,ux_avg_arr.imag,'-',color = 'b'); 
    plt.tight_layout();
    
    plt.subplot(325);
    for j in range(len(pars_1)): plt.plot(pars_1,Om_z_arr.real,'-',color = 'b'); 
#    plt.xlabel(key_2)
    plt.title('Om_z  ')
    plt.tight_layout();
    
    plt.subplot(326);
    for j in range(len(pars_1)): plt.plot(pars_1,Om_z_arr.imag,'-',color = 'b'); 
#    plt.xlabel(key_2)
    gridreso_SI = Header_Info['gridreso_SI']
    plt.title('$\Delta\Gamma$/n {:.1e}  {:.1e} gr {:.1f}'.format(pars_1[0],pars_1[-1],gridreso_SI*1e9 ))
    plt.tight_layout();
    plt.rcParams["figure.figsize"] = (6,4)
    plt.show()

    #coverage plots moved up
    if key_1 == "coverage":
        plt.subplot(131);
        plt.semilogx(pars_1,(dfc_arr.real + 716.8) / pars_1,'-')
        plt.xlabel(key_1); plt.title('$\Delta$f/n / coverage')
        
        plt.subplot(132);
        plt.semilogx(pars_1,(dfc_arr.imag - 716.8)  / pars_1,'-') 
        plt.xlabel(key_1); plt.title('$\Delta\Gamma$/n / coverage')
       
        plt.subplot(133);
        plt.plot(pars_1,((dfc_arr.imag-716.8)/(-dfc_arr.real-716.8)),'-') 
        plt.xlabel(key_1); plt.title('$\Delta\Gamma$/(-$\Delta$f)')
        plt.tight_layout();
        plt.show()


#    plt.subplot(131); plt.plot(pars_1,dfc_arr.real,'-'); 
#    plt.xlabel(key_1)
#    plt.title('df/n ' + key_1)
#    plt.subplot(132); plt.plot(pars_1,dfc_arr.imag,'-'); 
#    plt.xlabel(key_1)
#    plt.title('dG/n {:.1e}  {:.1e}'.format(pars_1[0],pars_1[-1]) )
#    plt.subplot(133); plt.plot(pars_1,comp_time_arr,'-'); 
#    plt.title('comp_time [s]')
#    plt.xlabel(key_1)
#    plt.show();
    data = [[0,0,0]]
    for i in range(len(dfc_arr)): 
        data = np.append(data,[[pars_1[i],dfc_arr.real[i],dfc_arr.imag[i]]], axis = 0)
    print(data);    
    np.savetxt('dfc_arr.txt',data)
  
    
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
    ux_avg_arr_12 = np.zeros((n1,n2),dtype = complex)
    Om_z_arr_12 = np.zeros((n1,n2),dtype = complex)

    for i_sim in range(n1 * n2):
        i1 = i_sim % n1;
        i2 = np.int(i_sim/n1);
        dfc_arr_12[i1,i2] = dfc_arr[i_sim]
        ux_avg_arr_12[i1,i2] = ux_avg_arr[i_sim]
        Om_z_arr_12[i1,i2] = Om_z_arr[i_sim]
    
    if key_2 == "cell_diam_by_sphere_diam":
        key_2 = "coverage"
        coverage = pars_2**-2 * np.pi / 4.; 
        pars_2 = coverage
   
    plt.subplot(231);
    for j in range(len(pars_1)): plt.plot(pars_2,dfc_arr_12[j,:].real,'-',color = Colors[j])
    plt.xlabel(key_2); plt.title('df/n ' + key_1)
    
    plt.subplot(234);
    for j in range(len(pars_1)): plt.plot(pars_2,dfc_arr_12[j,:].imag,'-',color = Colors[j]) 
    plt.xlabel(key_2); gridreso_SI = Header_Info['gridreso_SI']
    plt.title('$\Delta\Gamma$/n {:.1e}  {:.1e} gr {:.1f}'.format(pars_1[0],pars_1[-1],gridreso_SI*1e9 ))
    
    plt.subplot(232);
    for j in range(len(pars_1)): plt.plot(pars_2,ux_avg_arr_12[j,:].real - 1,'-',color = Colors[j]); 
    plt.xlabel(key_2); plt.title('ux_avg ' + key_1)
    
    plt.subplot(235);
    for j in range(len(pars_1)): plt.plot(pars_2,ux_avg_arr_12[j,:].imag,'-',color = Colors[j]);
    plt.xlabel(key_2)
    
    plt.subplot(233);
    for j in range(len(pars_1)): plt.plot(pars_2,Om_z_arr_12[j,:].real,'-',color = Colors[j]); 
    plt.xlabel(key_2)
    plt.title('Om_z  ' + key_1)
    
    plt.subplot(236);
    for j in range(len(pars_1)): plt.plot(pars_2,Om_z_arr_12[j,:].imag,'-',color = Colors[j]);
    plt.xlabel(key_2)
    
    plt.tight_layout();
    plt.show()

    if key_2 == "coverage":
        plt.subplot(221);
        for j in range(len(pars_1)): plt.plot(pars_2,(dfc_arr_12[j,:].real) / pars_2,'-',color = Colors[j])
        plt.xlabel(key_2); plt.title('($\Delta$f/n) / coverage')
        
        plt.subplot(222);
        for j in range(len(pars_1)): plt.plot(pars_2,(dfc_arr_12[j,:].imag)  / pars_2,'-',color = Colors[j]) 
        plt.xlabel(key_2); plt.title('($\Delta\Gamma$/n) / coverage')
       
        
        plt.subplot(223);
        for j in range(len(pars_1)): 
            plt.plot(pars_2,((dfc_arr_12[j,:].imag)/(-dfc_arr_12[j,:].real)),'-',color = Colors[j]) 
        plt.xlabel(key_2); plt.title('$\Delta\Gamma$/(-$\Delta$f) vs coverage')

        plt.subplot(224);
        for j in range(len(pars_1)): 
            plt.plot(-dfc_arr_12[j,:].real,((dfc_arr_12[j,:].imag)/(-dfc_arr_12[j,:].real)),'-',color = Colors[j]) 
        plt.xlabel('-$\Delta$f/n'); plt.title('$\Delta\Gamma$/(-$\Delta$f) vs -$\Delta$f')

        plt.tight_layout();
        plt.show()
        pars_2_even     = np.zeros(int(len(pars_2)/2.))
        pars_2_odd      = np.zeros(int(len(pars_2)/2.))
        dfc_arr_12_even = np.zeros((len(pars_1),int(len(pars_2)/2.)),dtype = complex)
        dfc_arr_12_odd  = np.zeros((len(pars_1),int(len(pars_2)/2.)),dtype = complex)
        for i in range(int(len(pars_2)/2.)):
            pars_2_even[i] = pars_2[2*i]
            pars_2_odd[ i] = pars_2[2*i + 1]
            for j in range(len(pars_1)):
                dfc_arr_12_even[j,i] = dfc_arr_12[j,2*i]
                dfc_arr_12_odd[ j,i] = dfc_arr_12[j,2*i+1]
        plt.subplot(221);
        for j in range(len(pars_1)): plt.plot(pars_2_even,(dfc_arr_12_even[j,:].real) / pars_2_even,'-',color = Colors[j])
        plt.xlabel(key_2); plt.title('($\Delta$f/n) / coverage even')
        
        plt.subplot(222);
        for j in range(len(pars_1)): plt.plot(pars_2_even,(dfc_arr_12_even[j,:].imag)  / pars_2_even,'-',color = Colors[j]) 
        plt.xlabel(key_2); plt.title('($\Delta\Gamma$/n) / coverage even')
        plt.subplot(223);
        for j in range(len(pars_1)): 
            plt.plot(pars_2_even,((dfc_arr_12_even[j,:].imag)/(-dfc_arr_12_even[j,:].real)),'-',color = Colors[j]) 
        plt.xlabel(key_2); plt.title('$\Delta\Gamma$/(-$\Delta$f) vs coverage')
        plt.subplot(224);
        for j in range(len(pars_1)): 
            plt.plot(-dfc_arr_12_even[j,:].real,((dfc_arr_12_even[j,:].imag)/(-dfc_arr_12_even[j,:].real)),'-',color = Colors[j]) 
        plt.xlabel('-$\Delta$f/n'); plt.title('$\Delta\Gamma$/(-$\Delta$f) vs -$\Delta$f')
        plt.tight_layout(); plt.show()
                
        plt.subplot(221);
        for j in range(len(pars_1)): plt.plot(pars_2_odd,(dfc_arr_12_odd[j,:].real) / pars_2_odd,'-',color = Colors[j])
        plt.xlabel(key_2); plt.title('($\Delta$f/n) / coverage odd')
        
        plt.subplot(222);
        for j in range(len(pars_1)): plt.plot(pars_2_odd,(dfc_arr_12_odd[j,:].imag)  / pars_2_odd,'-',color = Colors[j]) 
        plt.xlabel(key_2); plt.title('($\Delta\Gamma$/n) / coverage odd')
        plt.subplot(223);
        for j in range(len(pars_1)): 
            plt.plot(pars_2_odd,((dfc_arr_12_odd[j,:].imag)/(-dfc_arr_12_odd[j,:].real)),'-',color = Colors[j]) 
        plt.xlabel(key_2); plt.title('$\Delta\Gamma$/(-$\Delta$f) vs coverage')
        plt.subplot(224);
        for j in range(len(pars_1)): 
            plt.plot(-dfc_arr_12_odd[j,:].real,((dfc_arr_12_odd[j,:].imag)/(-dfc_arr_12_odd[j,:].real)),'-',color = Colors[j]) 
        plt.xlabel('-$\Delta$f/n'); plt.title('$\Delta\Gamma$/(-$\Delta$f) vs -$\Delta$f')
        plt.tight_layout(); plt.show()

    
    headerline = [0];
    for j in range(len(pars_1)) : headerline = np.append(headerline,pars_1[j]) 
    for j in range(len(pars_1)) : headerline = np.append(headerline,pars_1[j]) 
    print(headerline)
    data = [headerline]
    for i in range(len(pars_2)): 
        line = pars_2[i];
        for j in range(len(pars_1)) : line = np.append(line,dfc_arr_12[j,i].real) 
        for j in range(len(pars_1)) : line = np.append(line,dfc_arr_12[j,i].imag) 
        data = np.append(data,[line], axis = 0)
        print(line)
    np.savetxt('dfc_arr_12.txt',data)
        
    headerline = [0];
    for j in range(len(pars_1)) : headerline = np.append(headerline,pars_1[j]) 
    for j in range(len(pars_1)) : headerline = np.append(headerline,pars_1[j]) 
    print(headerline)
    data = [headerline]
    for i in range(len(pars_2_even)): 
        line = pars_2_even[i];
        for j in range(len(pars_1)) : line = np.append(line,dfc_arr_12_even[j,i].real) 
        for j in range(len(pars_1)) : line = np.append(line,dfc_arr_12_even[j,i].imag) 
        data = np.append(data,[line], axis = 0)
        print(line)
    np.savetxt('dfc_arr_12_even.txt',data)

    headerline = [0];
    for j in range(len(pars_1)) : headerline = np.append(headerline,pars_1[j]) 
    for j in range(len(pars_1)) : headerline = np.append(headerline,pars_1[j]) 
    print(headerline)
    data = [headerline]
    for i in range(len(pars_2_odd)): 
        line = pars_2_odd[i];
        for j in range(len(pars_1)) : line = np.append(line,dfc_arr_12_odd[j,i].real) 
        for j in range(len(pars_1)) : line = np.append(line,dfc_arr_12_odd[j,i].imag) 
        data = np.append(data,[line], axis = 0)
        print(line)
    np.savetxt('dfc_arr_12_odd.txt',data)
         
