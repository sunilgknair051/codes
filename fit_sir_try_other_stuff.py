import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = np.loadtxt('dfc_arr.txt', unpack=True)
df, dG = data[0], data[1] #, data[2]

i_discard = np.int(len(df)*2./3.)

time = time[i_discard:] - time[i_discard]
df = df[i_discard:]
dG = dG[i_discard:]

def func_exp(time, amp, tau, offset):
    return amp * np.exp(-time/tau) + offset

def func_dec_wave(time, amp, tau, offset,T,phi):
    return amp * np.exp(-time/tau)*np.cos(time/T + phi) + offset

#def func_
def func_parabola(time, amp, tau, offset):
    return amp*time**2 + tau*time + offset   

def Spiral(time, SpiralPars, offset):
#    offset = SpiralPars[0] + 1j*SpiralPars[1]
    ampl   = SpiralPars[2] + 1j*SpiralPars[3]
    omc    = SpiralPars[4] + 1j*SpiralPars[5]
    return ampl*np.exp(1j * omc * time) *offset

# Get parameters estimate
Guess_df = np.array([10,50])#,-460,35,0]) 
Guess_dG = np.array([10,50])#,1100,35,0]) 
#FitPars_df, pcov  = curve_fit(func_dec_wave, time, df, p0 = Guess_df)
#FitPars_dG, pcov = curve_fit(func_dec_wave, time, dG, p0 = Guess_dG)

FitPars_df, pcov  = curve_fit(Spiral, time, df)#, p0 = Guess_df)
FitPars_dG, pcov = curve_fit(Spiral,  time, dG)#, p0 = Guess_dG)

#popt
#df_fit = func_dec_wave(time,FitPars_df[0],FitPars_df[1],FitPars_df[2])#,FitPars_df[3],FitPars_df[4])
#dG_fit = func_dec_wave(time, FitPars_dG[0], FitPars_dG[1], FitPars_dG[2])#,FitPars_dG[3],FitPars_dG[4])

df_fit = Spiral(time,FitPars_df[0],FitPars_df[1],FitPars_df[2])#,FitPars_df[3],FitPars_df[4])
dG_fit = Spiral(time, FitPars_dG[0], FitPars_dG[1], FitPars_dG[2])#,FitPars_dG[3],FitPars_dG[4])

#print(test)
plt.plot(time,df,df_fit)
plt.show()
plt.plot(time,dG,'*',dG_fit)
plt.show()
print("Offset df {:.3f} Hz".format(FitPars_df[2]))
print("Offset dG {:.3f} Hz".format(FitPars_dG[2]))





