import numpy as np
import matplotlib.pyplot as plt
#          0    1    2    3    4    5     6
Colors = ['k', 'r', 'g', 'b', 'm', 'y', 'lime', \
 #          7                8        9          10          11        12
          'cornflowerblue','orange','indianred','chocolate','crimson','darkviolet'];
iovt = 2
i_ovt_str_arr = ['5 MHz', '15 MHz', '25 MHz', '35 MHz', '45 MHz', '55 MHz', '65 MHz'];

frac_above_surf_arr = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99,1.2,1.5,2.0]
dir_names = ['0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','0.99','1.2','1.5','2.0']
n_frac_above_surf = len(frac_above_surf_arr)
tandel_inv_arr = np.array([ 0.1     ,  0.158489,  0.251189,  0.398107,  0.630957,  1.      ,
        1.58489 ,  2.51189 ,  3.98107 ,  6.30957 , 10.      ])
nufac_arr = np.array([ 1.     ,   2.15443,   4.64159,  10.     ,  21.5443 ,  46.4159 ,
       100.     , 215.443  , 464.159])

dfcbyn_arr_12    = np.nan * np.ones((n_frac_above_surf,9,11),dtype = complex)
uxB_arr_12       = np.nan * np.ones((n_frac_above_surf,9,11),dtype = complex)
uxS_arr_12       = np.nan * np.ones((n_frac_above_surf,9,11),dtype = complex)
Rot_t_a_arr_12   = np.nan * np.ones((n_frac_above_surf,9,11),dtype = complex)
Shear_t_a_arr_12 = np.nan * np.ones((n_frac_above_surf,9,11),dtype = complex)
Dfcbyn_Part_arr_12  = np.nan * np.ones((n_frac_above_surf,9,11),dtype = complex)
Force_arr_12     = np.nan * np.ones((n_frac_above_surf,9,11),dtype = complex)
Torque_arr_12    = np.nan * np.ones((n_frac_above_surf,9,11),dtype = complex)
Txy_arr_12       = np.nan * np.ones((n_frac_above_surf,9,11),dtype = complex)
Dfbyn_Sph_arr_12 = np.nan * np.ones((n_frac_above_surf,9,11),dtype = complex)
Txy_Sph_W_by_CRad_arr_12 = np.nan * np.ones((n_frac_above_surf,9,11),dtype = complex)

uxB_arr_12[0,0] = np.nan + 1j * np.nan
Rot_t_a_arr_12[0,0] = np.nan + 1j * np.nan
Shear_t_a_arr_12[0,0] = np.nan + 1j * np.nan
Dfcbyn_Part_arr_12[0,0] = np.nan + 1j * np.nan

for i in range(n_frac_above_surf):
    data = np.loadtxt(dir_names[i] + '\\' + 'dfcbyn_21.txt',skiprows = 1);
    for j in range(11):
        dfcbyn_arr_12[i,:,j] = data[:,2*j+1] + 1j * data[:,2*j+2]   
    data = np.loadtxt(dir_names[i] + '\\' + 'uxB_21.txt',skiprows = 1);
    for j in range(11):
        uxB_arr_12[i,:,j] = data[:,2*j+1] + 1j * data[:,2*j+2]   
    data = np.loadtxt(dir_names[i] + '\\' + 'uxS_21.txt',skiprows = 1);
    for j in range(11):
        uxS_arr_12[i,:,j] = data[:,2*j+1] + 1j * data[:,2*j+2]   
    data = np.loadtxt(dir_names[i] + '\\' + 'Rot_t_a_21.txt',skiprows = 1);
    for j in range(11):
        Rot_t_a_arr_12[i,:,j] = data[:,2*j+1] + 1j * data[:,2*j+2]   
    data = np.loadtxt(dir_names[i] + '\\' + 'Shear_t_a_21.txt',skiprows = 1);
    for j in range(11):
        Shear_t_a_arr_12[i,:,j] = data[:,2*j+1] + 1j * data[:,2*j+2]   
    data = np.loadtxt(dir_names[i] + '\\' + 'Dfcbyn_Part_21.txt',skiprows = 1);
    for j in range(11):
        Dfcbyn_Part_arr_12[i,:,j] = data[:,2*j+1] + 1j * data[:,2*j+2]   
    data = np.loadtxt(dir_names[i] + '\\' + 'Force_21.txt',skiprows = 1);
    for j in range(11):
        Force_arr_12[i,:,j] = data[:,2*j+1] + 1j * data[:,2*j+2]   
    data = np.loadtxt(dir_names[i] + '\\' + 'Torque_21.txt',skiprows = 1);
    for j in range(11):
        Torque_arr_12[i,:,j] = data[:,2*j+1] + 1j * data[:,2*j+2]   
    data = np.loadtxt(dir_names[i] + '\\' + 'Txy_21.txt',skiprows = 1);
    for j in range(11):
        Txy_arr_12[i,:,j] = data[:,2*j+1] + 1j * data[:,2*j+2]   
    data = np.loadtxt(dir_names[i] + '\\' + 'Dfbyn_Sph_arr_12.txt',skiprows = 1);
    for j in range(11):
        Dfbyn_Sph_arr_12[i,:,j] = data[:,2*j+1] + 1j * data[:,2*j+2]   
    data = np.loadtxt(dir_names[i] + '\\' + 'Txy_Sph_W_by_CRad.txt',skiprows = 1);
    for j in range(11):
        Txy_Sph_W_by_CRad_arr_12[i,:,j] = data[:,2*j+1] + 1j * data[:,2*j+2]   

    uxB_arr_12[i,0,0] = np.nan + 1j * np.nan
    for j in range(5):
        for k in range(11):
            uxB_arr_12[-1,j,k] = uxB_arr_12[-2,j,k]
    Rot_t_a_arr_12[i,0,0] = np.nan + 1j * np.nan
    Shear_t_a_arr_12[i,0,0] = np.nan + 1j * np.nan
    Dfcbyn_Part_arr_12[i,0,0] = np.nan + 1j * np.nan


i_nufac = 2;
pars_plot = frac_above_surf_arr;
pars_1 = nufac_arr;
pars_2 = tandel_inv_arr;
#plt.figure(figsize=(8,9))
plt.rcParams["figure.figsize"] = (7,7)
plt.subplot(421);
for j in range(5): plt.plot(pars_plot,(dfcbyn_arr_12[:,i_nufac,2 * j].real),'.-',color = Colors[2*j])
plt.title('$\Delta$f/n) [Hz]');  ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(422);
for j in range(5): plt.plot(pars_plot,(dfcbyn_arr_12[:,i_nufac,2 * j].imag),'.-',color = Colors[2*j]) 
plt.title('$\Delta\Gamma$/n [Hz]'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(423);
for j in range(5): plt.plot(pars_plot,Dfbyn_Sph_arr_12[:,i_nufac,2 * j].real,'.-',color = Colors[2*j])
plt.title('$\Delta$f/n Sphere [Hz]'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(424);
for j in range(5): plt.plot(pars_plot,Dfbyn_Sph_arr_12[:,i_nufac,2 * j].imag,'.-',color = Colors[2*j])
plt.title('$\Delta\Gamma$/n Sphere [Hz]'); ax = plt.gca(); ax.axes.xaxis.set_visible(False); 
plt.subplot(425);
for j in range(5): plt.plot(pars_plot,Dfcbyn_Part_arr_12[:,i_nufac,2 * j].real,'.-',color = Colors[2*j])
plt.title('Df/n from Mass'); plt.xlabel('frac_above_surf');
plt.subplot(426);
for j in range(5): plt.plot(pars_plot,Dfcbyn_Part_arr_12[:,i_nufac,2 * j].imag,'.-',color = Colors[2*j],\
              label = str(tandel_inv_arr[2*j]))
plt.title('DG/n from Diss'); plt.xlabel('frac_above_surf');
plt.legend();

plt.subplot(427);
#label = str('$\eta$ [mPa s]' +  nufac_arr[i_nufac]) + 'width = 30, Rad = 5, height = 13'
label = '$\eta$ [mPa s] ' +  str(nufac_arr[i_nufac]) + '\n' + ' width = 30 nm \n Rad = 5 nm, \n height = 13 nm'\
       + '\n' + i_ovt_str_arr[iovt]\
       + '\n' + '(tan $\delta$)$^{-1}$ varied'
plt.plot([0,.1],[0,1],color = 'w',label = label); plt.legend();
ax = plt.gca(); ax.axes.xaxis.set_visible(False); ax.axes.yaxis.set_visible(False)
 


plt.savefig('vs_fac_above_surf1.png'); plt.tight_layout(); plt.show()


plt.subplot(421);
for j in range(5): plt.plot(pars_plot,np.abs(uxB_arr_12[:,i_nufac,2 * j]),'.-',color = Colors[2*j])
plt.title('uX Bottom Abs');  ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(422);
for j in range(5): plt.plot(pars_plot,np.angle(uxB_arr_12[:,i_nufac,2 * j])*180/np.pi,'.-',color = Colors[2*j]) 
plt.title('uX Bottom Phase'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(423);
for j in range(5): plt.plot(pars_plot,np.abs(uxS_arr_12[:,i_nufac,2 * j]),'.-',color = Colors[2*j])
plt.title('uX Avg'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(424);
for j in range(5): plt.plot(pars_plot,np.angle(uxS_arr_12[:,i_nufac,2 * j])*180/np.pi,'.-',color = Colors[2*j])
plt.title('uX Avg'); ax = plt.gca(); ax.axes.xaxis.set_visible(False); 
plt.subplot(425);
for j in range(5): plt.plot(pars_plot,np.abs(Rot_t_a_arr_12[:,i_nufac,2 * j]),'.-',color = Colors[2*j])
plt.title('Rotation x SphR'); ax = plt.gca(); ax.axes.xaxis.set_visible(False); 
plt.subplot(426);
for j in range(5): plt.plot(pars_plot,np.angle(Rot_t_a_arr_12[:,i_nufac,2 * j])*180/np.pi,'.-',color = Colors[2*j])
plt.title('Rotation x SphR'); ax = plt.gca(); ax.axes.xaxis.set_visible(False); 
plt.subplot(427);
for j in range(5): plt.plot(pars_plot,np.abs(Shear_t_a_arr_12[:,i_nufac,2 * j]),'.-',color = Colors[2*j])
plt.title('Shear x SphR'); plt.xlabel('frac_above_surf');
plt.subplot(428);
for j in range(5): plt.plot(pars_plot,np.angle(Shear_t_a_arr_12[:,i_nufac,2 * j])*180/np.pi,'.-',color = Colors[2*j])
plt.title('Shear x SphR'); plt.xlabel('frac_above_surf');

plt.savefig('vs_fac_above_surf2.png');
plt.tight_layout(); plt.show()

plt.subplot(421);
for j in range(5): plt.plot(pars_plot,np.abs(Force_arr_12[:,i_nufac,2 * j]),'.-',color = Colors[2*j])
plt.title('Force Abs');  ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(422);
for j in range(5): plt.plot(pars_plot,np.angle(Force_arr_12[:,i_nufac,2 * j])*180/np.pi,'.-',color = Colors[2*j]) 
plt.title('Force Phase'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(423);
for j in range(5): plt.plot(pars_plot,np.abs(Torque_arr_12[:,i_nufac,2 * j]),'.-',color = Colors[2*j])
plt.title('Torque Liq about Contact'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(424);
for j in range(5): plt.plot(pars_plot,np.angle(Torque_arr_12[:,i_nufac,2 * j])*180/np.pi,'.-',color = Colors[2*j])
plt.title('Torque Liq about Contact'); ax = plt.gca(); ax.axes.xaxis.set_visible(False); 
plt.subplot(425);
for j in range(5): plt.plot(pars_plot,np.abs(Txy_arr_12[:,i_nufac,2 * j]),'.-',color = Colors[2*j])
plt.title('Torque Liq about Cent'); ax = plt.gca(); ax.axes.xaxis.set_visible(False); 
plt.subplot(426);
for j in range(5): plt.plot(pars_plot,np.angle(Txy_arr_12[:,i_nufac,2 * j])*180/np.pi,'.-',color = Colors[2*j])
plt.title('Torque Liq about Cen'); ax = plt.gca(); ax.axes.xaxis.set_visible(False); 
plt.subplot(427);
for j in range(5): plt.plot(pars_plot,np.abs(Txy_Sph_W_by_CRad_arr_12[:,i_nufac,2 * j]),'.-',color = Colors[2*j])
plt.title('Torque @ Contact / ContRad'); plt.xlabel('frac_above_surf');
plt.subplot(428);
for j in range(5): plt.plot(pars_plot,np.angle(Txy_Sph_W_by_CRad_arr_12[:,i_nufac,2 * j])*180/np.pi,'.-',color = Colors[2*j])
plt.title('Torque @ Contact / ContRad'); plt.xlabel('frac_above_surf');

plt.savefig('vs_fac_above_surf3.png');
plt.tight_layout(); plt.show()

