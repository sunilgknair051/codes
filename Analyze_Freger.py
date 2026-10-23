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
	

def Read_Displacement_Sph(filename):
    data = np.transpose(np.loadtxt(filename))
    gridreso_SI = np.float(Header_Info['gridreso_SI']);
    rho_SI = 1e3;
    Mass_per_site = rho_SI * gridreso_SI**3;
    
    x = data[0];
    y = data[1];
    z = data[2];
    ux = data[3] + 1j * data[4]
    uy = data[5] + 1j * data[6]
    uz = data[7] + 1j * data[8]
    
    x_avg_CM = np.average(x);
    y_avg_CM = np.average(y); 
    z_avg_CM = np.average(z);
    
    print(x_avg_CM,y_avg_CM,z_avg_CM)
    Mom_Inertia_SI = Mass_per_site * np.sum((x-x_avg_CM)**2 + (y-y_avg_CM)**2) * gridreso_SI**2;
    print (Mom_Inertia_SI)
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
    return ux_avg, uy_avg, uz_avg, Om_x, Om_y, Om_z     


keys = Get_Labels()
Show_Plots = True;
header_filename = 'test083';
n_sims = 16;
HInfo_filename             = os.getcwd() + '\\' + 'results/' + header_filename +  '_HeaderInfo_' + str(n_sims).zfill(3) + '.dat';
Displacement_Sph_filename  = os.getcwd() + '\\' + 'results/' + header_filename + '_Displ_Sph_'   + str(n_sims).zfill(3) + '.dat';
Header_Info                   = Read_Header_Info(HInfo_filename,keys)
ux_avg, uy_avg, uz_avg, Om_x, Om_y, Om_z = Read_Displacement_Sph(Displacement_Sph_filename);
     
print(' ux = {:.3f} + i {:.3f} \n uy = {:.3f} + i {:.3f} \n uz = {:.3f} + i {:.3f}' \
      .format(ux_avg.real,ux_avg.imag,uy_avg.real,uy_avg.imag,uz_avg.real,uz_avg.imag))
print(' Om_x = {:.2e} + i {:.2e} \n Om_y = {:.2e} + i {:.2e} \n Om_z = {:.2e} + i {:.2e}' \
      .format(Om_x.real,Om_x.imag,Om_y.real,Om_y.imag,Om_z.real,Om_z.imag))
print(Header_Info)


#nufac = Header_Info['nufac']
#plt.subplot(222);
#for j in range(len(pars_1)):
#        plt.plot(ux,'-',color = 'b')#dfc_arr_12[:-1,     j].real,'-',color = 'b'); 
#        plt.xlabel(key_1)
#        plt.title(' ux vs.'+ key_1)#vs. nufac)#.format(pars_2[0],pars_2[-1]) )
##    plt.tight_layout();
#plt.show()
    #plt.title('dG/n {:.1e}  {:.1e} gr {:.1f}'.format(pars_1[0],pars_1[-1],gridreso_SI*1e9 ))
#plt.tight_layout();
#plt.show()
#plt.subplot(131); 
#plt.plot(Om_z,nufac,'-'); 
#plt.xlabel(Om_z)
#plt.title('ux_avg')
    #plt.subplot(132); plt.plot(); 
    #plt.xlabel()


