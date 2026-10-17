!Fortran code for particle near surface.
!S.G.

implicit none

! define constants 
real, parameter :: nu = 1/6  !0.166666
real, parameter :: pi = 3.1415927 
real, parameter :: Zq_SI = 8.8E+6
real, parameter :: nu_SI = 1E-6
real, parameter :: rho0 = 1
real, parameter :: rho_Liq_SI = 1000
real, parameter :: f0_SI = 5E+6
real, parameter :: delta_SI = 250E-9
real, parameter :: delfac = 1/6 !0.16666

real, parameter :: delta = NY*delfac
real, parameter :: omega = 2*nu / delta**2
real, parameter :: tauinv = 2/(6*nu+1)
real, parameter :: omtauinv = 1-tauinv
real, parameter :: uxWall = 1
!real, parameter :: 

!define weights and direction probabilities
 

subroutine init_equilibrium(dr,ux,uy)
!h = np.zeros((NX,NY,ndir),dtype = complex)
 for(int x = 0; x < Nx; ++x) 
	 {
       for(int y = 0; y < Ny ; ++y) 
	    {
	     for(int i = 0; i < ndir ; ++i) 
	      {
           cidotu = dirx[i]*ux[x,y] + diry[i]*uy[x,y]
           h[x,y,i] = wi[i] * dr[x,y]*(1 + 3*cidotu)
          }
        } 
	 }
	return h
end subroutine init_equilibrium

subroutine stream_h (h_src)
!h_dst = np.zeros((NX,NY,ndir),dtype = complex)
!FxLiqOnWall = 0; FxLiqOnSphe = 0;
  for(int x = 0; x < Nx; ++x) 
	 
       for(int y = 0; y < Ny ; ++y) 
	    
	     for(int i = 0; i < ndir ; ++i) 
	      
           IF inLiquid(x,y) THEN
        ibar = noslip[i]; xmd = (NX+x-dirx[i]) % NX; ymd = y-diry[i]
		     IF inLiquid(xmd, ymd) THEN
			     h_dst[x,y,i] = h_src[xmd,ymd,i]
		       ELSE IF  aboveTop(ymd) THEN
		         h_dst[x,y,i] = h_src[x,y,ibar]
		       ELSE IF belowBottom(ymd)
                     h_dst[x,y,i] = h_src[x,y,ibar] - 6 * wi[ibar]*dirx[ibar]
                 FxLiqOnWall += (h_dst[x,y,i] + h_src[x,y,ibar])*dirx[i]
			   ELSE inSample(xmd,ymd)
			     h_dst[x,y,i] = h_src[x,y,ibar] - 6 * wi[ibar] * dirx[ibar]
                 FxLiqOnSphe += (h_dst[x,y,i] + h_src[x,y,ibar])*dirx[i]
		     END IF
		   END IF	 
          
return h_dst,FxLiqOnWall,FxLiqOnSphe
end subroutine stream_h

subroutine compute_dr_u(h_in)
for(int x = 0; x < Nx; ++x) 
	 {
       for(int y = 0; y < Ny ; ++y) 
	    {
	       IF inLiquid(x,y) THEN
		     dr[x,y] = 0; ux[x,y] = 0; uy[x,y] = 0;
		    
			for(int i = 0; i < ndir ; ++i)
		       {
			     dr[x,y] += h_in[x,y,i]                   
                 ux[x,y] += h_in[x,y,i]*dirx[i]                   
                 uy[x,y] += h_in[x,y,i]*diry[i]
		 
                }
            END IF
	  	}	
	 }	
return dr,ux,uy
end subroutine compute_dr_u

subroutine collide(h_in,dr,ux,uy)
!h_out = np.ones((NX,NY,ndir),dtype = complex)
for(int x = 0; x < Nx; ++x) 
	 {
       for(int y = 0; y < Ny ; ++y) 
	    {
	       IF inLiquid(x,y) THEN
		     for(int i = 0; i < ndir ; ++i)
		      {
			    cidotu = dirx[i]*ux[x,y] + diry[i]*uy[x,y]
                    heq = wi[i]*(dr[x,y] + 3.*cidotu)   
                    h_out[x,y,i] = omtauinv * h_in[x,y,i] +
					tauinv * heq + 1j*omega * h_in[x,y,i] - omega**2/2 * heq # eq. 3.3
			  }
           END IF
		}
      }		
return h_out
end subroutine collide

subroutine calc_dfc()
TanStress = (FxLiqOnWall + FxLiqOnSphe)/ NX;
ZL = TanStress / uxWall; 
Zq = (rho0/rho_Liq_SI)*(nu/nu_SI)*(delta_SI/delta)*Zq_SI 
dfc = f0_SI* 1j /pi * ZL/Zq
dfc = -np.conjugate(dfc)
return dfc
end subroutine calc_dfc

end