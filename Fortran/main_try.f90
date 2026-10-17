!Couette Flow using LBM!

Program Couette_flow

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
!real,parameter :: uLowWallx = u_max*5
!real,parameter :: uLowWally = 0
!real,parameter :: uUpWallx = 0
!real,parameter :: uUpWally = 0
real,parameter  :: r = 1
real,parameter :: dummy,n,NSTEPS, h


subroutine main

f1 = init_equilibrium(rho0,ux,uy)
DO n = 1, NSTEPS
  dummy = stream(f1)
  print *, "dummy"
  (r,ux,uy) = compute_rho_u(f2)
  !print ()
  !print
  !print
  f1 = collide(f2,r,ux,uy)
!  IF n in (0,10,99) THEN
     !plot()
  ENDIF
ENDDO

end subroutine main

end Program Couette_flow
