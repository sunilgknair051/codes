import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from scipy.optimize import curve_fit


def taylor_green_inner(t,x,y):
    kx = 2. * np.pi/NX;
    ky = 2. * np.pi/NY;
    td = 1./(nu*(kx**2 + ky**2))
    X = x + 0.5
    Y = y + 0.5
    ux = - u_max * np.sqrt(ky/kx)*np.cos(kx*X)*np.sin(ky*Y)*np.exp(-1.*t/td)
    uy =   u_max * np.sqrt(kx/ky)*np.sin(kx*X)*np.cos(ky*Y)*np.exp(-1.*t/td)
    P = -0.25*rho0*u_max**2*((ky/kx)*np.cos(2.*kx*X) + (kx/ky)*np.cos(2.*ky*Y))*np.exp(-2.*t/td)
    ux = 0; uy = 0; P = 0;    
    rho = rho0 + 3. *P
    return rho,ux,uy      

def taylor_green_outer(t):
    for x in range(NX):
        for y in range(NY):
            r[x,y], ux[x,y], uy[x,y] = taylor_green_inner(t,x,y)
    return r, ux, uy

def arrows_plot(): 
    X, Y = np.meshgrid(np.arange(0, NX), np.arange(0, NY))
    plt.figure(1)
    plt.title('Taylor-Green')
    plt.quiver(Y, X, ux, uy)
    
def init_equilibrium(r,ux,uy):
    f = np.ones((NX,NY,ndir))
    for x in range(NX):
        for y in range(NY):
            for i in range(ndir):
                cidotu = dirx[i]*ux[x,y] + diry[i]*uy[x,y]
                if i_linearize == 0:
                    f[x,y,i] = wi[i]*r[x,y]*(1. + 3.*cidotu + 4.5*cidotu**2 - 1.5*(ux[x,y]**2+uy[x,y]**2))
                elif i_linearize == 1:     
                    f[x,y,i] = wi[i]*r[x,y]*(1. + 3.*cidotu) # run till x^1 to check decay
    return f

def insphere(x,y):
   for x in range(NX):
        for y in range(1,NY-1):
            if not insphere(x,y):
                    for i in range(ndir):
                        xmd = (NX+x-dirx[i]) % NX
                        ymd = (NY+y-diry[i]) % NY
                
                    if not insphere(xmd,ymd):
                         f_dst[x,y,i] = f_src[xmd,ymd,i]
                    else :
                         ibar = pn[i]
                         f_dst[x,y,i] = f_src[x,y,ibar] - 6 * wi[ibar]*(dirx[ibar]*Vel_bottom_X + diry[ibar]*Vel_bottom_Y)
                         
def stream(f_src):
    f_dst = np.ones((NX,NY,ndir))
    for x in range(NX):
        for y in range(1,NY-1):
            for i in range(ndir):
                xmd = (NX+x-dirx[i]) % NX
                ymd = (NY+y-diry[i]) % NY
                f_dst[x,y,i] = f_src[xmd,ymd,i]
 # 6  2  5      ^ y
# \ | /        |
# 3-0-1        | 
# / | \        | 
# 7  4  8       -----> x       
        
        y = 0
        #Bottom Wall Populations         #upward bounce
        for i in [2,5,6]:
            #cidotu = dirx[i]*Vel_bottom_X + diry[i]*Vel_bottom_Y   
            #print i
            ibar = pn[i]
            f_dst[x,y,i] = f_src[x,y,ibar] - 6 * wi[ibar]*(dirx[ibar]*Vel_bottom_X + diry[ibar]*Vel_bottom_Y)
            
        for i in [0,1,3,4,7,8]:              #for i in [0]and[1]....will not work!!
            xmd = (NX+x-dirx[i]) % NX
            ymd = (NY+y-diry[i]) % NY
            f_dst[x,y,i] = f_src[xmd,ymd,i]
#            print i
            
        y = NY - 1
        #Top Wall Populations           #downward bounce
        for i in [4,7,8]:
#            cidotu = dirx[i]*Vel_top_X + diry[i]*Vel_top_Y     
            ibar = pn[i]
            f_dst[x,y,i] = f_src[x,y,ibar] - 6 * wi[ibar]*(dirx[ibar]*Vel_top_X + diry[ibar]*Vel_top_Y)
            #print i
            
        for i in [0,1,2,3,5,6]:
            xmd = (NX+x-dirx[i]) % NX
            ymd = (NY+y-diry[i]) % NY
            f_dst[x,y,i] = f_src[xmd,ymd,i]
            #print i
        
         ### 
    

   
    
            

    return ((x-x_sphere)**2 + (y-y_sphere)**2 <= rad_sphere )
    return f_dst

def compute_rho_u(f):
    for x in range(NX):
        for y in range(NY):
            r[x,y] = 0;
            ux[x,y] = 0;
            uy[x,y] = 0;
            for i in range(ndir):
                r[x,y]  += f[x,y,i]                   
                ux[x,y] += f[x,y,i]*dirx[i]                   
                uy[x,y] += f[x,y,i]*diry[i]
            ux[x,y] = ux[x,y]/r[x,y]                    
            uy[x,y] = uy[x,y]/r[x,y]
    return r,ux,uy

def collide(f,r,ux,uy):
    f_out = np.ones((NX,NY,ndir))
    tauinv = 2./(6.*nu+1.)
    omtauinv = 1.-tauinv
    for x in range(NX):
        for y in range(NY):
            for i in range(ndir):
                cidotu = dirx[i]*ux[x,y] + diry[i]*uy[x,y]
                if i_linearize == 0:
                    feq = wi[i]*r[x,y]*(1. + 3.*cidotu + 4.5*cidotu**2 - 1.5*(ux[x,y]**2+uy[x,y]**2))
                elif i_linearize == 1:     
                    feq = wi[i]*r[x,y]*(1. + 3.*cidotu)
                f_out[x,y,i] = omtauinv * f[x,y,i] + tauinv * feq
    return f_out


# fit to an exponential     
def func(time, a, b):
    return a * np.exp(-time/b)

scale = 1; NX= 16*scale; NY = NX; t= 0; ndir = 9; w0 = 4./9.; ws = 1./9.; wd = 1./36.
i_linearize = 1
NSTEPS = 200* scale**2
wi = [w0,ws,ws,ws,ws,wd,wd,wd,wd]
dirx = [0,1,0,-1, 0,1,-1,-1, 1]
diry = [0,0,1, 0,-1,1, 1,-1,-1]
pn =   [0,3,4,1,2,7,8,5,6] 
# 6  2  5      ^ y
# \ | /        |
# 3-0-1        | 
# / | \        | 
# 7  4  8       -----> x
nu = 1./6.
tau = 3. * nu + 0.5
rho0 = 1;
u_max = 0.01/scale

#Declaring Wall Velocities 

Vel_top_X = 1.
Vel_top_Y = 0.
Vel_bottom_X = 0.
Vel_bottom_Y = 0.


#Sphere centre, velocity and radius values
x_sphere = 8.
y_sphere = 8.
vel_sphere = 5.
rad_sphere = 3.

dt = 0.

#circle_position = [8., 8., 2.]
##Defining a circle with an array
#def(circle_position):
#    xc = circles[0] #x-co of circle (center)
#    yc = circles[1] #y-co of circle (center)
#    r = circles[2] #radius of circle
#    arr=[]
#    for i in range(360):
#        y = yc + r*math.cos(i)
#        x = xc+ r*math.cos(i)
#        x=int(x)
#        y=int(y)
#        #Create array with all the x-co and y-co of the circle
#        arr.append([x,y])
#        circles=plt.Circle((x, y), r, color='b', fill=False)
#        plt.show()
#    return arr
##Initialized rho, ux and uy
#r = np.ones((NX,NY)); ux = np.zeros((NX,NY)); uy = np.zeros((NX,NY))
##r, ux, uy = taylor_green_outer(0)
#f1 = init_equilibrium(r,ux,uy)
#arrows_plot() 
#time = np.empty(0)
#u_decay = np.empty(0)

#f = open('filll.txt', 'wb')
for n in range(NSTEPS):
    f2 = stream(f1)
    r,ux,uy = compute_rho_u(f2)
#    time = np.append(time, n)
#    u_decay = np.append(u_decay, ux[10,15])
    #if np.mod(n,20) == 0:
#    print("%d %10.4e" % (n,ux[10,15])) #assign to an array
#    f.write("%10.4e %10.4e\n" % (n, ux[10,15]))
    f1 = collide(f2,r,ux,uy)
#f.close() 
arrows_plot() 
#plt.plot(xc, yc, '*')
circle=plt.Circle((x_sphere, y_sphere), rad_sphere, color='b', fill=False)
plt.axis('scaled')
fig = plt.gcf()
fig.gca().add_artist(circle)
#plt.plot(a_circle)
#plt.axis([0, 16, 0, 16])
#plt.subplot(111)
#plt.plot(xc,yc, '*')   
plt.show()

fig.savefig('venn_diagramm.png')
#data = np.loadtxt('filll.txt', unpack=True)
#time, ydata = data[0], data[1]

"""
popt, pcov = curve_fit(func, time, ydata)

test = func(time,popt[0],popt[1])

#print(test)
plt.figure(2)
plt.plot (test)
plt.plot(time,ydata)
"""
