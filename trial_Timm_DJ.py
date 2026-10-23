# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 08:45:32 2018

@author: User
"""
import numpy as np


#LBM in python D2Q9#

#c = np.array([(x,y) for x in [0,-1,1] for y in [0,-1,1]]) # Lattice velocities.
#print (c)
#w_i = weights...define
#c_iu = c_i dot u   

#Variables involved#



 

#Initialization of particle population with equilibrium values#
def equilibrium (f, r, u, v):
     for y in xrange(0, NY):
        for x in xrange(0, NX):
            rho = r
            ux  = u 
            uy  = v 
            
            for i in xrange(0, ndir):
                c_iu = x *ux + y * uy
                f_i = w_i*rho*(1 + 3*c_iu + 4.5*(c_iu)^2 -1.5*(ux^2-uy^2))
                
        

#Particle Streaming Step#

def stream (f_src, f_dst):
    for y in xrange(0, NY):
        for x in xrange(0, NX):
            for i in xrange(0, ndir):
                
                xmd = (NX+x-dirx)%NX
                ymd = (NY+y-diry)%NY
                
                f_dist  = f_src;

#Computation of density and velocity#
                
def compute_rho (f, r, u, v):
    for y in xrange(0, NY):
        for x in xrange(0, NX):
            rho = 0
            ux  = 0
            uy  = 0
            
        for i in xrange(0, ndir):
            rho  = rho + f
            ux   = ux + dirx*f
            uy   = uy + diry*f
         
            r = rho
            u = ux/rho
            v = uy/rho
                
#Particle Collision Step#
            
def collide (f, r, u, v):
 
    tau_inv = 2/(6*nu+1)
    m_tau_inv = 1 - tau_inv
    
    for y in xrange(0, NY):
        for x in xrange(0, NX):
            rho = r
            ux = u
            uy = v
            
            
            for i in xrange(0, ndir):
                c_iu = dirx*ux + diry*uy
                f_eq = w_i*rho*(1 + 3*c_iu + 4.5*(c_iu)^2 -1.5*(ux^2-uy^2))
                f = m_tau_inv *f + tau_inv*f_eq
                
                
#Main function call#

def main (f1 , f2):
  
   equilibrium (f1, rho,  ux, uy)

for n in xrange(0, NSTEPS):
    stream (f1, f2)
    compute_rho (f2, rho, ux, uy)
    collide (f2, rho, )                       
                
                
            








