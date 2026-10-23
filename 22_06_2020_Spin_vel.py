from __future__ import division
import xdrlib
import math
import numpy as np
import glob
import matplotlib.pyplot as plt
from matplotlib import cm 
import os.path
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

    f = open('results/test072_Displ_Sph_058.dat','r')
    lines = f.readlines()
    columns = lines[0].split()
#print (lines)
    x = int(columns[0])
    y = int(columns[1])
    z = int(columns[2])
    ux = np.zeros((x,y,z),dtype = complex)
    uy = np.zeros((x,y,z),dtype = complex)
    uz = np.zeros((x,y,z),dtype = complex)

#print (lines)
#print (len(lines))
    gridreso_SI = np.float(Header_Info['gridreso_SI']);
    rho_SI = 1e3;
    Mass_per_site = rho_SI * gridreso_SI**3;

    for i in range(len(lines)):
        columns = lines[i].split()
#   x = int(columns[0]) 
#   y = int(columns[1]) 
#   z = int(columns[2])
#   print (x,y,z)
#   print(columns[0],columns[1],columns[2],columns[3],\
#         columns[4],columns[5],columns[6],columns[7],columns[8])
        ux = float(columns[3]) + 1j * float(columns[4])
#   print(ux)
        uy = float(columns[5]) + 1j * float(columns[6])
        uz = float(columns[7]) + 1j * float(columns[8])
   
        x_avg_CM = np.average(x);
#   print(x_avg_CM)
        y_avg_CM = np.average(y); 
        z_avg_CM = np.average(z);
        Mom_Inertia_SI = Mass_per_site * np.sum((x-x_avg_CM)**2 + (y-y_avg_CM)**2) * gridreso_SI**2;
#   print (str(Mom_Inertia_SI))  need to print the whole 10e-39 value
        ux_avg = np.average(ux);
        uy_avg = np.average(uy);
        uz_avg = np.average(uz); 
   
        L_x =( Mass_per_site * (uy * (z - z_avg_CM) - uz * (y - y_avg_CM)) * gridreso_SI**2)
        print(L_x)
        L_y =( Mass_per_site * (uz * (x - x_avg_CM) - ux * (z - z_avg_CM)) * gridreso_SI**2)
        L_z =( Mass_per_site * (ux * (y - y_avg_CM) - uy * (x - x_avg_CM)) * gridreso_SI**2)
    # Rotation
        Om_x = np.sum((L_x)) / Mom_Inertia_SI 
#   print(Om_x)
        Om_y = np.sum((L_y)) / Mom_Inertia_SI 
        Om_z = np.sum((L_z)) / Mom_Inertia_SI 
        f.close
        return ux_avg, uy_avg, uz_avg, Om_x, Om_y, Om_z 
    
keys = Get_Labels()
Show_Plots = True;
header_filename = 'test072';
n_sims = 76;
HInfo_filename             = os.getcwd() + '\\' + 'results/' + header_filename +  '_HeaderInfo_' + str(n_sims).zfill(3) + '.dat';
Displacement_Sph_filename  = os.getcwd() + '\\' + 'results/' + header_filename + '_Displ_Sph_'   + str(n_sims).zfill(3) + '.dat';
Header_Info                   = Read_Header_Info(HInfo_filename,keys)
ux_avg, uy_avg, uz_avg, Om_x, Om_y, Om_z = Read_Displacement_Sph(Displacement_Sph_filename);



nufac = Header_Info['nufac']
plt.subplot(222);
for j in range(len(pars_1)):
        plt.plot(ux,'-',color = 'b')#dfc_arr_12[:-1,     j].real,'-',color = 'b'); 
        plt.xlabel(key_1)
        plt.title(' ux vs.'+ key_1)#vs. nufac)#.format(pars_2[0],pars_2[-1]) )
        plt.tight_layout();
plt.show()
    #plt.title('dG/n {:.1e}  {:.1e} gr {:.1f}'.format(pars_1[0],pars_1[-1],gridreso_SI*1e9 ))
#plt.tight_layout();
#plt.show()
#plt.subplot(131); 
#plt.plot(Om_z,nufac,'-'); 
#plt.xlabel(Om_z)
#plt.title('ux_avg')
    #plt.subplot(132); plt.plot(); 
    #plt.xlabel()























     
#print(' ux = {:.3f} + i {:.3f} \n uy = {:.3f} + i {:.3f} \n uz = {:.3f} + i {:.3f}' \
#      .format(ux_avg.real,ux_avg.imag,uy_avg.real,uy_avg.imag,uz_avg.real,uz_avg.imag))
#print(' Om_x = {:.2e} + i {:.2e} \n Om_y = {:.2e} + i {:.2e} \n Om_z = {:.2e} + i {:.2e}' \
#      .format(Om_x.real,Om_x.imag,Om_y.real,Om_y.imag,Om_z.real,Om_z.imag))
#print(Header_Info)

 




#   ux[x,y,z] = float(columns[3]) + 1j * float(columns[4])
#   uy[x-1,y-1,z-1] = float(columns[5]) + 1j * float(columns[6])
#   uz[x-1,y-1,z-1] = float(columns[7]) + 1j * float(columns[8])
#   ux[i-1] = float(columns[0])
#   uy[i-1] = float(columns[1])
#   uz[i-1] = float(columns[2])
   #print (x,y,z)
#   print (ux,uy,uz)

#f = plt.figure(figsize=(15,15))
#dropinifac = 1/5.
#plt.subplot(221); plt.plot(-dfs[int(dropinifac*len(dfs)):], '-'); plt.plot(dGs[int(dropinifac*len(dfs)):], '-');
#plt.subplot(222); plt.plot(-dfs[int(dropinifac*len(dfs)):],dGs[int(dropinifac*len(dfs)):], '-'); 
##f = plt.figure(figsize=(15,15))
#plt.subplot(223); plt.imshow(ux.transpose().real,cmap=cm.Reds);
#plt.subplot(224); plt.imshow(uy.transpose().real,cmap=cm.Reds);

#plt.legend()
#plt.title("Freqs")
#plt.xlabel("T")
#plt.ylabel("F")
#plt.savefig("foo" + ".pdf", format='pdf')
#plt.show()

#
#plt.subplot(); plt.plot(y_arr,ux.real,'-'); plt.plot(y_arr,ux.imag,'-'); 

#plt.subplot(421); plt.imshow(uxsecz.real.transpose(),cmap=cm.Reds); 
#plt.subplot(422); plt.imshow(uysecz.real.transpose(),cmap=cm.Reds);
#plt.subplot(423); plt.imshow(uxsecy.real.transpose(),cmap=cm.Reds);
#plt.subplot(424); plt.imshow(uysecy.real.transpose(),cmap=cm.Reds);
#plt.subplot(425); plt.plot(ux[:,int(ny/2.0),int(nz/2.0)].real,'-'); 
#plt.subplot(426); plt.plot(ux[:,int(ny/2.0),int(nz/2.0)].imag,'-'); 
#plt.subplot(427); plt.plot(ux[int(nx/2.0),int(ny/2.0),:].real,'-'); 
#plt.subplot(428); plt.plot(ux[int(nx/2.0),int(nz/2.0),:].imag,'-'); 
    



	



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






