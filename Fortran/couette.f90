!Couette Flow using LBM!

implicit none

!Declare constants!

real,parameter :: scale = 1
real,parameter :: NX = 8*scale
real,parameter :: NY = NX
real,parameter :: t = 0
real,parameter :: ndir = 9
real,parameter :: w0 = 4/9
real,parameter :: ws = 1/9
real,parameter :: wd = 1/36
real,parameter :: nu = 1/6
real,parameter :: tau = 3*nu + 0.5
real,parameter :: rho0 = 1
real,parameter :: NSTEPS = 100*scale**2
real,parameter :: u_max = 0.04/scale
real,parameter :: uLowWallx = u_max*5
real,parameter :: uLowWally = 0
real,parameter :: uUpWallx = 0
real,parameter :: uUpWally = 0







subroutine init_equilibrium (r,ux,uy)
! f = np.ones ((NX, NY, ndir))

for (int x=0; x <Nx ; ++x)
  for (int y=0; y<Ny; ++y)
      for (int i=0; i<ndir; ++i)
          cidotu = dirx(i) * ux(x,y) + diry(i) * uy(x,y)
          f(x,y,i) = wi(i)*r(x,y)(1+3*cidotu + 4.5*cidotu**2  - 1.5 *(ux(x,y)**2 + uy(x,y)**2))

return f
end subroutine init_equilibrium



!subroutine stream_pop (f_src)
!f_dst = np.ones((NX,NY,ndir))
for (int x=0; x <Nx ; ++x)
  for (int y=1; y<Ny; ++y)
      for (int i=0; i<ndir; ++i)
         xmd = (NX+x-dirx(i))%NX
         ymd = (NY+y-diry(i))%NY
         f_dst(x,y,i)  = f_src(xmd,ymd,i)
      end







          cidotu = dirx(i) * ux(x,y) + diry(i) * uy(x,y)
          f(x,y,i) = wi(i)*r(x,y)(1+3*cidotu + 4.5*cidotu**2  - 1.5 *(ux(x,y)**2 + uy(x,y)**2))



return f_dst
end subroutine stream_pop



subroutine compute_rho_u
for (int x=0; x <Nx ; ++x)
  for (int y=0; y<Ny; ++y)
      for (int i=0; i<ndir; ++i)
          cidotu = dirx(i) * ux(x,y) + diry(i) * uy(x,y)
          f(x,y,i) = wi(i)*r(x,y)(1+3*cidotu + 4.5*cidotu**2  - 1.5 *(ux(x,y)**2 + uy(x,y)**2))


return r,ux,uy
end subroutine compute_rho_u



subroutine collide(f,r,ux.uy)
for (int x=0; x <Nx ; ++x)
  for (int y=0; y<Ny; ++y)
      for (int i=0; i<ndir; ++i)
          cidotu = dirx(i) * ux(x,y) + diry(i) * uy(x,y)
          f(x,y,i) = wi(i)*r(x,y)(1+3*cidotu + 4.5*cidotu**2  - 1.5 *(ux(x,y)**2 + uy(x,y)**2))
return f_out
end subroutine collide


end
