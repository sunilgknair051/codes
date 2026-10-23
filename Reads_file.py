from __future__ import division
#import xdrlib
#import math
import numpy as np
import glob
import matplotlib.pyplot as plt



files = glob.glob('rep1text.asc')
files.sort()
##print len(files)
mddata = np.zeros((len(files),10))

for ti, filen in enumerate(files):
	f = open(filen,'r')
	c = 0
	for l in f:
		data = [float(a) for a in l.split()]
		mddata[ti] = data
		c += 1

nx = mddata[:,0]
nx = np.trim_zeros(nx,'b') 
ny = mddata[:,1]
ny = np.trim_zeros(ny,'b') 

print (nx, ny)

plt.plot(nx,ny,label='force x-direction',marker='o')
#plt.plot(time,Stokes,label='Stokes Drag',marker='o')
#
#plt.legend()
#plt.title("D")
#plt.xlabel("T")
#plt.ylabel("F")
#plt.savefig("h" + ".pdf", format='pdf')
#plt.show()
#
#
