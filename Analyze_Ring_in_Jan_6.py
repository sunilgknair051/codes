import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import os
import configparser
Show_Plots = True;
#          0    1    2    3    4    5     6
Colors = ['k', 'r', 'g', 'b', 'm', 'y', 'lime', \
 #          7                8        9          10          11        12
          'cornflowerblue','orange','indianred','chocolate','crimson','darkviolet'];
i_ovt_str_arr = ['5 MHz', '15 MHz', '25 MHz', '35 MHz', '45 MHz', '55 MHz', '65 MHz'];
def Read_Header_Info(filename):
    keys =  ["ncycles_max", "rad_Sph_SI", "cellheight_SI", \
	"gridreso_SI", 	"cell_diam_by_sphere_diam", "frac_above_surf", \
	"nufac", "tandelinv", "iovt", "nx", "ny", "nz", "ncycles_max", "outInterval","jqmax",\
    "SphR","Sphx","Sphy","Sphz","t_computation_secs"];

    config = configparser.ConfigParser(); config.read(filename); sectionname = header_filename;
    Header_Info = dict(zip(keys,np.zeros(len(keys))))
    for key in keys:                                                                              
        if key in ['ncycles_max','iovt', 'nx', 'ny', 'nz', 'outInterval']:
            Header_Info[key] = int(config.get(sectionname, key,fallback=-1))
        else :
            Header_Info[key] = float(config.get(sectionname, key,fallback=-9999.))
    nx = np.int(Header_Info['nx']);        
    ny = np.int(Header_Info['ny']);        
    nz = np.int(Header_Info['nz']);        
    iovt = np.int(Header_Info['iovt']);    
    SphR = np.float(Header_Info['SphR']);
    Sphx = np.float(Header_Info['Sphx']);
    Sphy = np.float(Header_Info['Sphy']);
    Sphz = np.float(Header_Info['Sphz']);
    nufac =  np.float(Header_Info['nufac'])        
    tandelinv  = np.float(Header_Info['tandelinv'])
    frac_above_surf = np.float(Header_Info['frac_above_surf'])            
    comp_time = np.float(Header_Info['t_computation_secs'])
    return keys,Header_Info,nx,ny,nz,iovt,SphR,Sphx,Sphy,Sphz,nufac,tandelinv,frac_above_surf,comp_time 

def Analyze_Pops(filename,dfbyn_extrapol,dGbyn_extrapol) :
    h = np.zeros((nx,ny,nz,ndir),dtype = complex);
    data = np.loadtxt(filename); 
    for j in range(len(data)): 
        x = int(data[j,0]); y = int(data[j,1]); z = int(data[j,2]); i = int(data[j,3]);
        h[x,y,z,i] = data[j,4] + 1j * data[j,5];
    dr  = np.zeros((nx,ny,nz),dtype = complex)
    ux  = np.zeros((nx,ny,nz),dtype = complex)
    uy  = np.zeros((nx,ny,nz),dtype = complex)
    uz  = np.zeros((nx,ny,nz),dtype = complex)
    for x in range (nx):
        for y in range (ny):
            for z in range (nz):
                for i in range (ndir):
                    dr[x,y,z]  += h[x,y,z,i];
                    ux[x,y,z]  += dirx[i] * h[x,y,z,i];
                    uy[x,y,z]  += diry[i] * h[x,y,z,i];
                    uz[x,y,z]  += dirz[i] * h[x,y,z,i];
    xCM = 0; yCM = 0; zCM = 0; norm = 0; ACont = 0;
    uxS = 0; uyS = 0; uzS = 0; uxB = 0;
    for x in range (nx):
        for y in range (ny):
            for z in range (nz):
                r = ((x - Sphx)**2 + (y - Sphy)**2 + (z - Sphz)**2)**0.5;
                if r < SphR : 
                    xCM += x; yCM += y; zCM += z;  
                    uxS += ux[x,y,z]; uyS += uy[x,y,z]; uzS += uz[x,y,z]; 
                    norm += 1; 
                    if y == 0 : 
                        ACont += 1;
                        uxB += ux[x,y,z]; 
    if norm > 0:
        Volume = norm;
        xCM /= norm; yCM /= norm; zCM /= norm;  
        uxS /= norm; uyS /= norm; uzS /= norm;
    rCont = (ACont / np.pi)**0.5; 
    if ACont > 0 : uxB /= ACont    
    eijS = np.zeros((3,3), dtype = complex); norm = 0;    
    for x in range (nx):
        for y in range (ny):
            for z in range (nz):
                r = ((x - Sphx)**2 + (y - Sphy)**2 + (z - Sphz)**2)**0.5;
                if r < SphR : 
                    xm = (nx+x-1) % nx; 
                    xp = (nx+x+1) % nx; 
                    eijS[0,0] += (ux[xp,y,z] - ux[xm,y,z])/2. 
                    eijS[1,0] += (uy[xp,y,z] - uy[xm,y,z])/2. 
                    eijS[2,0] += (uz[xp,y,z] - uz[xm,y,z])/2. 
                    norm += 1;

                    ym = y-1
                    yp = y+1
                    if ym >= 0 : 
                        eijS[0,1] += (ux[x,yp,z] - ux[x,ym,z])/2. 
                        eijS[1,1] += (uy[x,yp,z] - uy[x,ym,z])/2. 
                        eijS[2,1] += (uz[x,yp,z] - uz[x,ym,z])/2.
                    else :    
                        eijS[0,1] += (ux[x,yp,z] - 1)/2. 
                        eijS[1,1] += (uy[x,yp,z] - 0)/2. 
                        eijS[2,1] += (uz[x,yp,z] - 0)/2.
                    norm += 1;

                    zm = (nz+z-1) % nz; 
                    zp = (nz+z+1) % nz; 
                    eijS[0,2] += (ux[x,y,zp] - ux[x,y,zm])/2. 
                    eijS[1,2] += (uy[x,y,zp] - uy[x,y,zm])/2. 
                    eijS[2,2] += (uz[x,y,zp] - uz[x,y,zm])/2. 
                    norm += 1;
    if norm > 0 : eijS /= norm;
#    Omx = eijS[1,2] - eijS[2,1] 
#    Omy = eijS[2,0] - eijS[0,2] 
    Omz = 1./2.*(eijS[0,1] - eijS[1,0]) 
#    Syz = eijS[1,2] + eijS[2,1] 
#    Szx = eijS[2,0] + eijS[0,2] 
    Sxy = 1./2.*(eijS[0,1] + eijS[1,0]) 
#    Trace = eijS[0,0] + eijS[1,1] + eijS[2,2]
#    Sxx = eijS[0,0] - Trace / 3. 
#    Syy = eijS[1,1] - Trace / 3. 
#    Com = Trace;
    Rot_t_a = Omz * yCM;
    Shear_t_a = Sxy * yCM;
    
    TijS = np.zeros((3,3), dtype = complex); norm = 0;
    Fx = 0; Fy = 0; Fz = 0; Fx_Sph_On_Wall = 0; Txy_Sph_On_Wall = 0; 
   
                
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                r = ((x - Sphx)**2 + (y - Sphy)**2 + (z - Sphz)**2)**0.5;
                if r < SphR : 
                    for i in range(ndir):
                        ibar = noslip[i];
                        xmd = (nx+x-dirx[i]) % nx; 
                        ymd = y-diry[i]
                        zmd = (nz+z-dirz[i]) % nz; 
                        rmd = ((xmd - Sphx)**2 + (ymd - Sphy)**2 + (zmd - Sphz)**2)**0.5;
                        if rmd >= SphR : 
                            x2 = 1./2.*(x + xmd) 
                            y2 = 1./2.*(y + ymd) 
                            z2 = 1./2.*(z + zmd)
                            Fx +=  (h[x,y,z,i] + h[x,y,z,ibar]) * dirx[i]
                            Fy +=  (h[x,y,z,i] + h[x,y,z,ibar]) * diry[i]
                            Fz +=  (h[x,y,z,i] + h[x,y,z,ibar]) * dirz[i]
                            TijS[0,0] += h[x,y,z,i] * dirx[i] * (x2 - xCM)
                            TijS[1,0] += h[x,y,z,i] * diry[i] * (x2 - xCM)
                            TijS[2,0] += h[x,y,z,i] * dirz[i] * (x2 - xCM)
                            TijS[0,1] += h[x,y,z,i] * dirx[i] * (y2 - yCM)
                            TijS[1,1] += h[x,y,z,i] * diry[i] * (y2 - yCM)
                            TijS[2,1] += h[x,y,z,i] * dirz[i] * (y2 - yCM)
                            TijS[0,2] += h[x,y,z,i] * dirx[i] * (z2 - zCM)
                            TijS[1,2] += h[x,y,z,i] * diry[i] * (z2 - zCM)
                            TijS[2,2] += h[x,y,z,i] * dirz[i] * (z2 - zCM)
                            norm += 1;
                        if ymd < 0 :
                            h2 = h[x,y,z,ibar] - 6.0 * wi[ibar] * dirx[ibar];    
                            Fx_Sph_On_Wall  += (h2 + h[x,y,z,ibar]) * dirx[i];
                            Txy_Sph_On_Wall += (h2 + h[x,y,z,ibar]) * diry[i] * (x2 - Sphx);

    if norm > 0 : TijS /= norm; 
#    Tx = TijS[1,2] - TijS[2,1] 
#    Ty = TijS[2,0] - TijS[0,2] 
#    Tz = TijS[0,1] - TijS[1,0] 
#    Tyz = TijS[1,2] + TijS[2,1] 
#    Tzx = TijS[2,0] + TijS[0,2] 
    Txy = TijS[0,1] + TijS[1,0] 
#    Trace = TijS[0,0] + TijS[1,1] + TijS[2,2]
#    Txx = TijS[0,0] - Trace / 3. 
#    Tyy = TijS[1,1] - Trace / 3. 
#    TCom = Trace;
    if frac_above_surf < 1 :
        Txy_Sph_On_Wall_by_CRad = Txy_Sph_On_Wall / (SphR * (1 - (frac_above_surf-0.5)**2)**0.5 );
    else:
        Txy_Sph_On_Wall_by_CRad = np.nan;
    Force = Fx;
    Torque = Txy + Fx * yCM;
    f0_SI = 5e6; dq_SI = 332.075e-6; 
    ovtorder = 2 * iovt + 1; nu = 1./6.
    rholiq = 1; rhoq = 2.65; 
    gridreso_SI = np.float(Header_Info['gridreso_SI']);
    nu_SI = 1e-6; delta_SI = (2.0 * nu_SI / (2.0 *  np.pi * f0_SI))**0.5;
    delta = delta_SI / gridreso_SI;
    omega = 2 * nu / delta**2
    f0 = omega / (2 * np.pi * ovtorder);

    dq = dq_SI / gridreso_SI
#    Zq = f0 * rhoq * dq / 2.  Here was the mistake
    Zq = 2. * f0 * rhoq * dq
    Dfbyn_Mass =  - omega * rholiq * Volume / (nx * nz)/ (np.pi * Zq) * f0_SI / ovtorder; 
    if dfbyn_extrapol < -0.5 :
        TrapLiqFrac = (dfbyn_extrapol - Dfbyn_Mass) / dfbyn_extrapol  + 1j * 0.
    else :    
        TrapLiqFrac = 0  + 1j * 0.
    if tandelinv > 0: 
        tandel = 1 / tandelinv
#        cos = 1.     / (1 + tandel**2)**0.5;
        sin = tandel / (1 + tandel**2)**0.5;
        nu_Sph_re = nu * nufac * sin;
    else : nu_Sph_re = nu * nufac
    Qdot = 0;
    for i in range(3):
        for j in range(3):
            Qdot += nu_Sph_re * Volume * np.abs((eijS[i,j] + eijS[j,i])/2.)**2
    ERes = Zq / (4. * f0)
#    print('Qdot',Qdot,'Zq',Zq,'f0',f0,'eijS[0,1]' ,eijS[0,1] )     
    DGbyn_Diss = Qdot / (nx * nz) / ERes /(4. * np.pi) * f0_SI / f0 / ovtorder
#    print('DG_Diss',DG_Diss)     
    uxB = uxB - 1;
    uxS = uxS - 1;
    Dfcbyn_Part = Dfbyn_Mass + 1j * DGbyn_Diss
    Dfbyn_Sph = f0_SI * 1j / (np.pi* Zq) * Fx_Sph_On_Wall / (nx * nz) / ovtorder;
    Deriv_Pars = np.array([uxB,uxS,Rot_t_a,Shear_t_a,Dfcbyn_Part,Force,Torque,Txy,\
                           Dfbyn_Sph,Txy_Sph_On_Wall_by_CRad,TrapLiqFrac])
    return dr,ux,uy,uz,xCM,yCM,zCM,rCont,Deriv_Pars
	
def Write_fit_Results():
    fname = os.getcwd() + '\\' + 'results/' + header_filename + '_dfcs.txt'
    f = open(fname,'w+'); f.close();
    f = open(fname,'a+'); 
    for i in range(len(pars_1)):
        for j in range(len(pars_2)):
            print('{:d} {:d} {:.2e} {:.2e} {:.2f}  {:.2f}'.\
              format(i, j, pars_1[i], pars_2[j], \
                     dfcbyn_arr_12[i,j].real, dfcbyn_arr_12[i,j]).real);
            f.write('{:d} {:d} {:.2e} {:.2e} {:.2f} \t {:.2f} \n'.\
              format(i, j, pars_1[i], pars_2[j], \
                     dfcbyn_arr_12[i,j].real, dfcbyn_arr_12[i,j].imag))      
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

def Fit_Spiral(filename):
    data = np.loadtxt(filename); 
    dfbyns = data[:,0]; dGbyns = data[:,1]; t_arr = np.arange(len(dfbyns));
    discardini = int(len(dfbyns)*1./3.)
    t_arr = t_arr[discardini:] - t_arr[discardini] 
    dfbyns = dfbyns[discardini:] 
    dGbyns = dGbyns[discardini:]
    dfcbyn = dfbyns + 1j * dGbyns
    if len(dfbyns) > 8 :
        SpiralParsGuess = np.array([dfbyns[-1],dGbyns[-1], dfbyns[0]-dfbyns[-1],dGbyns[0]-dGbyns[-1],\
                                    0.1/len(t_arr),2./len(t_arr)]) 
    
        SpiralParsFit,dfcbyn_fitted,chi2 = fit_model(t_arr, dfcbyn,SpiralParsGuess)
        dfbyn_extrapol = SpiralParsFit[0]
        dGbyn_extrapol = SpiralParsFit[1]
        residual = dfcbyn - dfcbyn_fitted
        if Show_Plots and (n_sims % gap_display == 0) and (ref == "no") : 
            plt.rcParams["figure.figsize"] = (7,2)
            plt.subplot(221); plt.plot(t_arr,-dfbyns,'-'); plt.plot(t_arr,dGbyns,'-'); 
            plt.plot(t_arr,-dfcbyn_fitted.real,'-'); plt.plot(t_arr,dfcbyn_fitted.imag,'-');
            plt.title(("df/n:%7.1f  dG/n:%7.1f "%(dfbyn_extrapol,dGbyn_extrapol)))
            plt.subplot(222); plt.plot(residual.real,'-'); plt.plot(residual.imag,'-'); 
            plt.title(("residuals" + str(n_sims)))
            plt.tight_layout(); plt.show()
            if n_sims == 21 : 
                data = np.array([dfbyns,dGbyns,dfcbyn_fitted.real,dfcbyn_fitted.imag])
                np.savetxt('spiral21.txt',np.transpose(data));
            if n_sims == 29 : 
                data = np.array([dfbyns,dGbyns,dfcbyn_fitted.real,dfcbyn_fitted.imag])
                np.savetxt('spiral29.txt',np.transpose(data));
    else : dfbyn_extrapol,dGbyn_extrapol,chi2 = np.nan,np.nan,np.nan       
    return dfbyn_extrapol,dGbyn_extrapol,chi2

def write_to_file(fname,data_array):
    headerline = [0];
    for j in range(len(pars_1)) : 
        headerline = np.append(headerline,pars_1[j]) 
        headerline = np.append(headerline,pars_1[j]) 
    data = [headerline]
    for i in range(len(pars_2)): 
        line = pars_2[i];
        for j in range(len(pars_1)) : 
            line = np.append(line,data_array[j,i].real) 
            line = np.append(line,data_array[j,i].imag) 
        data = np.append(data,[line], axis = 0)
    np.savetxt(fname + '.txt',data)

ndir = 19; zW = 1./3.; aW = 1./18.; dW = 1./36.;
#         0 1  2 3  4 5  6 7  8 9  10 11 12 13 14 15 16 17 18
dirx   = [0,1,-1,0, 0,0, 0,1,-1,1 ,-1, 0, 0, 1,-1, 1,-1, 0, 0];
diry   = [0,0,0 ,1,-1,0, 0,1,-1,0 , 0, 1,-1,-1, 1, 0, 0, 1,-1];
dirz   = [0,0,0 ,0, 0,1,-1,0, 0,1 ,-1, 1,-1, 0, 0,-1, 1,-1, 1];
noslip = [0,2,1 ,4, 3,6, 5,8, 7,10, 9,12,11,14,13,16,15,18,17];
wi = [zW, aW,aW,aW,aW,aW,aW, dW,dW,dW,dW,dW,dW,dW,dW,dW,dW,dW,dW];

# header_filename = open('header_filename.txt', 'r').read();
header_filename = 'test005';
n_sims_max = 60;
gap_display = 1;
 
dfcbyn_arr = []; comp_time_arr = []; chi2_arr = []; 
Header_Info_arr = []; Deriv_Pars_arr = [np.zeros(11,dtype = complex)]
uxS_arr = []; Om_z_arr = [];
n_sims = 0; i_file_exists = 1;
#while ((i_file_exists == 1) and n_sims < 1000):
while ((i_file_exists == 1) and n_sims < n_sims_max):
    if n_sims % 10 == 0: print(n_sims);
#    HInfo_filename        = os.getcwd() + '\\' + 'results/' + header_filename + '_HeaderInfo_' + str(n_sims).zfill(3) + '.dat';
#    dfc_by_n_filename     = os.getcwd() + '\\' + 'results/' + header_filename + '_ring_in_'    + str(n_sims).zfill(3) + '.dat';
#    dfc_by_n_ref_filename = os.getcwd() + '\\' + 'results/' + header_filename + '_ring_in_'    + str(n_sims).zfill(3) + '_ref.dat';
#    Pop_filename          = os.getcwd() + '\\' + 'results/' + header_filename + '_Pop_'        + str(n_sims).zfill(3) + '.dat';
    HInfo_filename        = os.getcwd() + '\\' + 'results/' + header_filename + '_HeaderInfo_' + str(n_sims).zfill(3) + '.dat';
    dfc_by_n_filename     = os.getcwd() + '\\' + 'results/' + header_filename + '_ring_in_'    + str(n_sims).zfill(3) + '.dat';
    dfc_by_n_ref_filename = os.getcwd() + '\\' + 'results/' + header_filename + '_ring_in_'    + str(n_sims).zfill(3) + '_ref.dat';
    Pop_filename          = os.getcwd() + '\\' + 'results/' + header_filename + '_Pop_'        + str(n_sims).zfill(3) + '.dat';

    if os.path.isfile(HInfo_filename) :
        ref = "no"
        keys,Header_Info,nx,ny,nz,iovt,SphR,Sphx,Sphy,Sphz,nufac,tandelinv,frac_above_surf,comp_time = Read_Header_Info(HInfo_filename)
        dfbyn_extrapol,dGbyn_extrapol,chi2  = Fit_Spiral(dfc_by_n_filename);
        if os.path.isfile(dfc_by_n_ref_filename) : 
            ref = "yes"
            dfbyn_extrapol_ref,dGbyn_extrapol_ref,chi2  = Fit_Spiral(dfc_by_n_ref_filename);
        else : dfbyn_extrapol_ref,dGbyn_extrapol_ref = dfbyn_extrapol,dGbyn_extrapol    
        dfbyn_extrapol -= dfbyn_extrapol_ref;
        dGbyn_extrapol -= dGbyn_extrapol_ref;
        if os.path.isfile(Pop_filename) :
            dr,ux,uy,uz,xCM,yCM,zCM,rCont,Deriv_Pars = Analyze_Pops(Pop_filename,dfbyn_extrapol,dGbyn_extrapol)
        
        Header_Info_arr = np.append(Header_Info_arr, [Header_Info],                    axis=0);  
        dfcbyn_arr      = np.append(dfcbyn_arr,      [dfbyn_extrapol + 1j * dGbyn_extrapol], axis=0);  
        chi2_arr        = np.append(chi2_arr,        [chi2],                           axis=0); 
        comp_time_arr   = np.append(comp_time_arr,   [comp_time],                      axis=0); 
        Deriv_Pars_arr  = np.append(Deriv_Pars_arr,  [Deriv_Pars],                     axis=0); 
        n_sims += 1;

#        if Show_Plots : 
#            x1 = np.arange(len(ux[0,:])); y1 = np.arange(len(ux[:,0])); X, Y = np.meshgrid(x1, y1)
#            plt.subplot(221); cp = plt.contourf(Y,X,ux.real); plt.colorbar(cp); 
#            plt.subplot(222); cp = plt.contourf(Y,X,ux.imag); plt.colorbar(cp); 
#            plt.subplot(223); cp = plt.contourf(Y,X,uy.real); plt.colorbar(cp); 
#            plt.subplot(224); cp = plt.contourf(Y,X,uy.imag); plt.colorbar(cp); 
#            plt.tight_layout(); plt.show();
    else : 
        i_file_exists = 0;
        print('file not found: ' + HInfo_filename)

key_1 = 'none'; key_2 = 'none';
if n_sims > 1 : 
    for key in keys:
        if Header_Info_arr[1][key] != Header_Info_arr[0][key] and \
          key != 't_computation_secs' and \
          key != 'cell_diam_by_sphere_diam' and key != 'nz' and key != 'Sphx'  and key != 'Sphz': 
              key_1 = key

for key in keys:
    pars = np.zeros(n_sims)
    for i_sim in range(n_sims) : pars[i_sim] = Header_Info_arr[i_sim][key]
    if (np.std(pars) > 1e-12) and key != key_1 and key != 't_computation_secs' and \
    key != 'cell_diam_by_sphere_diam' and key != 'nz' and key != 'Sphx'  and key != 'Sphz': 
        key_2 = key;
        
if key_2 == 'none': 
    n1 = n_sims;
    n2 = 1;
else :     
    y_start = Header_Info_arr[0][key_2]
    i_sim = 0;
    while i_sim < n_sims-1 : 
        i_sim += 1;
        diff = np.abs(Header_Info_arr[i_sim][key_1] - Header_Info_arr[0][key_1])
        if diff < 1e-12 : break
    n1 = i_sim;
    n2 = np.int((n_sims+1)/n1);

pars_1 = np.zeros(n1); pars_2 = np.zeros(n2)
for i1 in range(n1) : pars_1[i1] = Header_Info_arr[i1]    [key_1]
if key_2 != 'none' : 
    for i2 in range(n2) : pars_2[i2] = Header_Info_arr[i2*n1] [key_2]
else : pars_2[0] = np.nan     
dfcbyn_arr_12 = np.zeros((n1,n2),dtype = complex)
Deriv_Pars_arr_12 = np.zeros((n1,n2,11),dtype = complex)

for i_sim in range(n1 * n2):
    i1 = i_sim % n1;
    i2 = np.int(i_sim/n1);
    dfcbyn_arr_12[i1,i2] = dfcbyn_arr[i_sim]
    Deriv_Pars_arr_12[i1,i2] = Deriv_Pars_arr[i_sim] 

uxB_arr_12       = Deriv_Pars_arr_12[:,:,0]
uxS_arr_12       = Deriv_Pars_arr_12[:,:,1]
Rot_t_a_arr_12   = Deriv_Pars_arr_12[:,:,2]
Shear_t_a_arr_12 = Deriv_Pars_arr_12[:,:,3]
Dfcbyn_Part_arr_12  = Deriv_Pars_arr_12[:,:,4]
Force_arr_12     = Deriv_Pars_arr_12[:,:,5]
Torque_arr_12    = Deriv_Pars_arr_12[:,:,6]
Txy_arr_12       = Deriv_Pars_arr_12[:,:,7]
Dfbyn_Sph_arr_12 = Deriv_Pars_arr_12[:,:,8]
Txy_Sph_W_by_CRad_arr_12 = Deriv_Pars_arr_12[:,:,9]
TrapLiqFrac_12   = Deriv_Pars_arr_12[:,:,10]

#uxB_arr_12[0,0] = np.nan + 1j * np.nan
#Rot_t_a_arr_12[0,0] = np.nan + 1j * np.nan
#Shear_t_a_arr_12[0,0] = np.nan + 1j * np.nan
#Dfcbyn_Part_arr_12[0,0] = np.nan + 1j * np.nan

pars_2B = np.pi * 5**2 / pars_2**2 * 100
key_2B = 'coverage [%]'
   
plt.rcParams["figure.figsize"] = (7,9)
plt.subplot(421);
for j in range(len(pars_1)): plt.plot(pars_2B,(dfcbyn_arr_12[j,:].real),'-',color = Colors[j])
plt.title('$\Delta$f/n [Hz]');  ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(422);
for j in range(len(pars_1)): plt.plot(pars_2B,(dfcbyn_arr_12[j,:].imag),'-',color = Colors[j]) 
plt.title('$\Delta\Gamma$/n [Hz]'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(423);
for j in range(len(pars_1)): plt.plot(pars_2B,Dfbyn_Sph_arr_12[j,:].real,'-',color = Colors[j])
plt.title('$\Delta$f/n Sphere [Hz]'); plt.xlabel(key_2B);
plt.subplot(424);
for j in range(len(pars_1)): plt.plot(pars_2B,Dfbyn_Sph_arr_12[j,:].imag,'-',color = Colors[j])
plt.title('$\Delta\Gamma$/n Sphere [Hz]'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(425);
for j in range(len(pars_1)): plt.plot(pars_2B,TrapLiqFrac_12[j,:].real,'.-',color = Colors[j])
plt.title('Trapped Liq Frac'); plt.xlabel(key_2B);
plt.subplot(426);
for j in range(len(pars_1)): plt.plot(pars_2B,Dfcbyn_Part_arr_12[j,:].imag,'.-',color = Colors[j],\
              label = str(round(pars_1[j],3)))
plt.title('DG/n from Diss'); plt.xlabel(key_2B);
plt.legend();

plt.subplot(427);
cwdstr = os.getcwd(); cwdstr2 = cwdstr[-5:]     
#label = str('$\eta$ [mPa s]' +  nufac_arr[i_nufac]) + 'width = 30, Rad = 5, height = 13'
label = cwdstr2 + '\n' + ' width = 30 nm \n Rad = 5 nm, \n height = 13 nm'\
       + '\n' + i_ovt_str_arr[iovt]\
       + '\n' + '|$\eta$| varied'
plt.plot([0,.1],[0,1],color = 'w',label = label); plt.legend();
ax = plt.gca(); ax.axes.xaxis.set_visible(False); ax.axes.yaxis.set_visible(False)

plt.tight_layout(); plt.savefig('vs_nufac1.png'); plt.show()
 


plt.subplot(421);
for j in range(len(pars_1)): plt.plot(pars_2B,np.abs(uxB_arr_12[j,:]),'-',color = Colors[j])
plt.title('uX Bottom Abs'); plt.xlabel(key_2B);
plt.subplot(422);
for j in range(len(pars_1)): plt.plot(pars_2B,np.angle(uxB_arr_12[j,:])*180/np.pi,'-',color = Colors[j])
plt.title('uX Bottom Phase'); plt.xlabel(key_2B);
plt.subplot(423);
for j in range(len(pars_1)): plt.plot(pars_2B,np.abs(uxS_arr_12[j,:]),'-',color = Colors[j])
plt.title('uX Avg'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(424);
for j in range(len(pars_1)): plt.plot(pars_2B,np.angle(uxS_arr_12[j,:])*180/np.pi,'-',color = Colors[j])
plt.title('uX Avg'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(425);
for j in range(len(pars_1)): plt.plot(pars_2B,np.abs(Rot_t_a_arr_12[j,:]),'-',color = Colors[j])
plt.title('Rotation x SphR'); 
plt.subplot(426);
for j in range(len(pars_1)): plt.plot(pars_2B,np.angle(Rot_t_a_arr_12[j,:])*180/np.pi,'-',color = Colors[j])
plt.title('Rotation x SphR'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(427);
for j in range(len(pars_1)): plt.plot(pars_2B,np.abs(Shear_t_a_arr_12[j,:]),'-',color = Colors[j])
plt.title('Shear x SphR');  ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(428);
for j in range(len(pars_1)): plt.plot(pars_2B,np.angle(Shear_t_a_arr_12[j,:])*180/np.pi,'-',color = Colors[j])
plt.title('Shear x SphR'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.tight_layout(); plt.savefig('vs_nufac2.png'); plt.show()
     

plt.subplot(421);
for j in range(len(pars_1)): plt.plot(pars_2B,np.abs(Force_arr_12[j,:]),'-',color = Colors[j])
plt.title('Force Abs'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(422);
for j in range(len(pars_1)): plt.plot(pars_2B,np.angle(Force_arr_12[j,:])*180/np.pi,'-',color = Colors[j])
plt.title('Force Phase'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(423);
for j in range(len(pars_1)): plt.plot(pars_2B,np.abs(Torque_arr_12[j,:]),'-',color = Colors[j])
plt.title('Torque Liq about Contact');  ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(424);
for j in range(len(pars_1)): plt.plot(pars_2B,np.angle(Torque_arr_12[j,:])*180/np.pi,'-',color = Colors[j])
plt.title('Torque Liq about Contact'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(425);
for j in range(len(pars_1)): plt.plot(pars_2B,np.abs(Txy_arr_12[j,:]),'-',color = Colors[j])
plt.title('Torque Liq about Cent'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(426);
for j in range(len(pars_1)): plt.plot(pars_2B,np.angle(Txy_arr_12[j,:])*180/np.pi,'-',color = Colors[j])
plt.title('Torque Liq about Cent'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(427);
for j in range(len(pars_1)): plt.plot(pars_2B,np.abs(Txy_Sph_W_by_CRad_arr_12[j,:]),'-',color = Colors[j])
plt.title('Torque @ Contact / Crad'); plt.xlabel(key_2B); 
plt.subplot(428);
for j in range(len(pars_1)): plt.plot(pars_2B,np.angle(Txy_Sph_W_by_CRad_arr_12[j,:])*180/np.pi,'-',color = Colors[j])
plt.title('Torque @ Contact / CRad'); plt.xlabel(key_2B); 
plt.tight_layout(); plt.savefig('vs_nufac3.png'); plt.show()

plt.subplot(421);
for j in range(len(pars_2)): plt.plot(pars_1[:],(dfcbyn_arr_12[:,j].real),'-',color = Colors[j])
plt.title('$\Delta$f/n [Hz]');  ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(422);
for j in range(len(pars_2)): plt.plot(pars_1[:],(dfcbyn_arr_12[:,j].imag),'-',color = Colors[j]) 
plt.title('$\Delta\Gamma$/n [Hz]'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(423);
for j in range(len(pars_2)): plt.plot(pars_1[:],Dfbyn_Sph_arr_12[:,j].real,'-',color = Colors[j])
plt.title('$\Delta$f/n Sphere [Hz]'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(424);
for j in range(len(pars_2)): plt.plot(pars_1[:],Dfbyn_Sph_arr_12[:,j].imag,'-',color = Colors[j])
plt.title('$\Delta\Gamma$/n Sphere [Hz]'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(425);
for j in range(len(pars_2)): plt.plot(pars_1[:],TrapLiqFrac_12[:,j].real,'-',color = Colors[j])
plt.title('Trapped Liq Frac'); plt.xlabel(key_1); 
plt.subplot(426);
for j in range(len(pars_2)): plt.plot(pars_1[:],Dfcbyn_Part_arr_12[:,j].imag,'-',color = Colors[j],\
              label = str(round(pars_2[j],3)))
plt.title('DG/n from Diss'); plt.xlabel(key_1); 
plt.legend();

plt.subplot(427);
cwdstr = os.getcwd(); cwdstr2 = cwdstr[-5:]     
#label = str('$\eta$ [mPa s]' +  nufac_arr[i_nufac]) + 'width = 30, Rad = 5, height = 13'
label = cwdstr2 + '\n' + ' width = 30 nm \n Rad = 5 nm, \n height = 13 nm'\
       + '\n' + i_ovt_str_arr[iovt]\
       + '\n' + '(tan $\delta$)$^{-1}$ varied'
plt.plot([0,.1],[0,1],color = 'w',label = label); plt.legend();
ax = plt.gca(); ax.axes.xaxis.set_visible(False); ax.axes.yaxis.set_visible(False)
plt.tight_layout(); plt.savefig('vs_tandel1.png'); plt.show()


plt.subplot(421);
for j in range(len(pars_2)): plt.plot(pars_1[:],np.abs(uxB_arr_12[:,j]),'-',color = Colors[j])
plt.title('uX Bottom Abs'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(422);
for j in range(len(pars_2)): plt.plot(pars_1[:],np.angle(uxB_arr_12[:,j])*180/np.pi,'-',color = Colors[j])
plt.title('uX Bottom Phase'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(423);
for j in range(len(pars_2)): plt.plot(pars_1[:],np.abs(uxS_arr_12[:,j]),'-',color = Colors[j])
plt.title('uX Avg'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(424);
for j in range(len(pars_2)): plt.plot(pars_1[:],np.angle(uxS_arr_12[:,j])*180/np.pi,'-',color = Colors[j])
plt.title('uX Avg'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(425);
for j in range(len(pars_2)): plt.plot(pars_1[:],np.abs(Rot_t_a_arr_12[:,j]),'-',color = Colors[j])
plt.title('Rotation x SphR'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(426);
for j in range(len(pars_2)): plt.plot(pars_1[:],np.angle(Rot_t_a_arr_12[:,j])*180/np.pi,'-',color = Colors[j])
plt.title('Rotation x SphR'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(427);
for j in range(len(pars_2)): plt.plot(pars_1[:],np.abs(Shear_t_a_arr_12[:,j]),'-',color = Colors[j])
plt.title('Shear x SphR');  plt.xlabel(key_1); 
plt.subplot(428);
for j in range(len(pars_2)): plt.plot(pars_1[:],np.angle(Shear_t_a_arr_12[:,j])*180/np.pi,'-',color = Colors[j])
plt.title('Shear x SphR'); plt.xlabel(key_1); 
plt.tight_layout(); plt.savefig('vs_tandel2.png'); plt.show()

plt.subplot(421);
for j in range(len(pars_2)): plt.plot(pars_1[:],np.abs(Force_arr_12[:,j]),'-',color = Colors[j])
plt.title('Force Abs'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(422);
for j in range(len(pars_2)): plt.plot(pars_1[:],np.angle(Force_arr_12[:,j])*180/np.pi,'-',color = Colors[j])
plt.title('Force Phase'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(423);
for j in range(len(pars_2)): plt.plot(pars_1[:],np.abs(Torque_arr_12[:,j]),'-',color = Colors[j])
plt.title('Torque Liq about Contact');  ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(424);
for j in range(len(pars_2)): plt.plot(pars_1[:],np.angle(Torque_arr_12[:,j])*180/np.pi,'-',color = Colors[j])
plt.title('Torque Liq about Contact'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(425);
for j in range(len(pars_2)): plt.plot(pars_1[:],np.abs(Txy_arr_12[:,j]),'-',color = Colors[j])
plt.title('Torque Liq about Cent'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(426);
for j in range(len(pars_2)): plt.plot(pars_1[:],np.angle(Txy_arr_12[:,j])*180/np.pi,'-',color = Colors[j])
plt.title('Torque Liq about Cent'); ax = plt.gca(); ax.axes.xaxis.set_visible(False);
plt.subplot(427);
for j in range(len(pars_2)): plt.plot(pars_1[:],np.abs(Txy_Sph_W_by_CRad_arr_12[:,j]),'-',color = Colors[j])
plt.title('Torque @ Contact / CRad'); plt.xlabel(key_1); 
plt.subplot(428);
for j in range(len(pars_2)): plt.plot(pars_1[:],np.angle(Txy_Sph_W_by_CRad_arr_12[:,j])*180/np.pi,'-',color = Colors[j])
plt.title('Torque @ Contact / CRad'); plt.xlabel(key_1); 
plt.tight_layout(); plt.savefig('vs_tandel3.png')

plt.show()


write_to_file('dfcbyn_21',dfcbyn_arr_12)
write_to_file('uxB_21',uxB_arr_12)
write_to_file('uxS_21',uxS_arr_12)
write_to_file('Rot_t_a_21',Rot_t_a_arr_12)
write_to_file('Shear_t_a_21',Shear_t_a_arr_12)
write_to_file('Dfcbyn_Part_21',Dfcbyn_Part_arr_12)
write_to_file('TrapLiqFrac',TrapLiqFrac_12)
write_to_file('Force_21',Force_arr_12)
write_to_file('Torque_21',Torque_arr_12)
write_to_file('Txy_21',Txy_arr_12)
write_to_file('Dfbyn_Sph_arr_12',Dfbyn_Sph_arr_12)
write_to_file('Txy_Sph_W_by_CRad',Txy_Sph_W_by_CRad_arr_12)

