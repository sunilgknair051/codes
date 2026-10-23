import xdrlib
import math
import numpy as np
import glob
import matplotlib.pyplot as plt

#Velocity = 0.0015
#shear_rate = (0.0015/48)
#u_x = shear_rate * 0.75 
#volumefrac = ((4./3)*math.pi*125)/(48*48*48)
#Stokes_Drag = 6 * math.pi * (1/6)* 5 * u_x
#Hasimoto_Drag = 6 *math.pi*(1/6)* 5 *0.0015*(1+1.76(volumefrac**(1./3.)))
#print Stokes_Drag
#print Hasimoto_Drag


#files = glob.glob('output/md-cfg_default_t*-*.asc')
#files.sort()
##print len(files)
#mddata = np.zeros((len(files),100,23))
##print files
#for ti, filen in enumerate(files):
	#f = open(filen,'r')
	#c = 0
	#for lindex, l in enumerate(f):
		#data = [float(a) for a in l.split()]
		#mddata[ti,lindex] = data
		#c += 1
##print mddata[0]
##for file in files:
#xs = mddata[0,:,1] 
#xs = np.trim_zeros(xs,'b')

files = glob.glob('output/md-cfg_default_t*-*.asc')
files.sort()
print len(files)
mddata = np.zeros((len(files),23))
for ti, filen in enumerate(files):
	f = open(filen,'r')
	c = 0
	for l in f:
		data = [float(a) for a in l.split()]
		mddata[ti] = data
		c += 1

xs = mddata[:,3]
time = mddata[:,0]
#ys = mddata[:,8]
#xpos = [x*(48.0/50.0) for x in range(0,51)]
#print ys
#print xpos

#print len(xs)
#print len(xpos)

#xs = mddata[:,0,7]
#xpos = mddata[:,0,1]
#print xs;
#print mddata;
#np.average()
#ys = mddata[,8]
#zs = mddata[:,0,9]
#vel= mddata[:,0,4]
#Formula for Stokes Drag--> Fd=6(pi).eta.R.V
#F_d = ((math.pi)*5*2*0.0015*0.75)
#ar = np.ones(100)
#sh=F_d
#F_d = mddata[:,0,2.3e-2]
#print mddata
#print(xs)
#x = [1,100]0
#y = [0.023]

plt.plot(time,xs,label='xs',marker='o')
#plt.plot(xpos,xs,label='xs',marker='o')
#plt.plot(sh,label='stokes', marker='o')
#plt.plot(ys,label='ys',marker='o')
#plt.plot(zs,label='zs',marker='o')
#plt.plot(x, y)
plt.legend()
plt.title("YO!")
plt.xlabel("Time")
plt.ylabel("ANgVelocity")
plt.show()
plt.savefig("yes.pdf")

#import xdrlib
#import numpy as np
#import glob
#import matplotlib.pyplot as plt


#files = glob.glob('output/md-cfg_output_t*-*.asc')
#files.sort()

#mddata = np.zeros((len(files),64,20))

#for ti, filen in enumerate(files):
	#f = open(filen,'r')
	#c = 0
	#for l in f:
		#data = [float(a) for a in l.split()]
		#mddata[ti,c] = data
		#c += 1

#xs = mddata[:,0,1]
#ys = mddata[:,0,2]
#zs = mddata[:,0,3]

#print(xs)

#plt.plot(xs,label='xs',marker='o')
##plt.plot(ys,label='ys',marker='o')
##plt.plot(zs,label='zs',marker='o')
#plt.legend()
#plt.show()


