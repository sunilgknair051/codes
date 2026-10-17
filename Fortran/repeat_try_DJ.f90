!Couette Flow LBM!
!!The LBM code with D2Q9 solves for a couette flow with a u_max of 0.4!!


!this module declares the variables being used and we can just call the module throughout to use the variables!
module declare
    implicit none
    real :: w0,ws,wd
    parameter(w0=4./9., ws=1./9., wd= 1./36.)
    integer :: nx,ny,ndir
    parameter(nx=20,ny=20,ndir=8)
!    real :: f(1:nx,1:ny,0:ndir),fn(1:nx,1:ny,0:ndir),f_out(1:nx,1:ny,0:ndir)
    real :: f1main(1:nx,1:ny,0:ndir),f2main(1:nx,1:ny,0:ndir)
    real :: rmain(1:nx,1:ny),uxmain(1:nx,1:ny),uymain(1:nx,1:ny)
    integer :: dirx(0:8)=(/0,1,0,-1,0,1,-1,-1,1/)
    integer :: diry(0:8)=(/0,0,1,0,-1,1,1,-1,-1/)
    real    :: wi(0:8)=(/w0,ws,ws,ws,ws,wd,wd,wd,wd/)
    real,parameter :: nu = 1./6.
    real,parameter :: nu = 1./6.
    real,parameter :: tau = 3. * nu + 0.5
    real,parameter :: rho0 = 1.
    real,parameter :: tauinv = 2./(6.*nu + 1.)
    real,parameter :: omtauinv = 1. - tauinv
    real,parameter :: scale = 1. , u_max = 0.4/scale
    integer, parameter :: NSTEPS = 1000
    real,parameter :: uLowWallx = u_max*5., uLowWally = 0., uUpWallx = 0., uUpWally = 0.
    integer,parameter :: pn(0:8)=(/0,3,4,1,2,7,8,5,6/)
    integer,parameter :: out_unit = 20
end module declare

!The main program we call repeat as well as call all the functions being used!
program repeat
use declare
implicit none
integer :: x,y,i
DO x = 1,nx
  DO y = 1,ny
    rmain(x,y) = 1
    uxmain(x,y) = 0
    uymain(x,y) = 0
  ENDDO
ENDDO

call init_equilibrium(f1main,rmain,uxmain,uymain)
do i = 1,NSTEPS
    call stream(rmain,f1main,f2main)
    call compute_rho_u(f2main,rmain,uxmain,uymain)
    call collide(f2main,rmain,uxmain,uymain,f1main)
end do

!creating an ASCII file named rep1text.asc and saving details of the simulation!!
open(6,file='rep1text.asc')!, form = 'formatted')
!write(6,*)          nx,ny
write(6,*)'          x     ','     y    ','ux(x,y)','         uy(x,y)'
DO x= 1,nx
  DO y =1,ny
    write(6,*)h2main
  ENDDO
ENDDO
close(6)
end program repeat

!!The LBM involves steps such as equilibrium, streaming, updating/computing Rho and, collisions!!
!!An important point to observe here is the equation we use is only upto the order 2.!!
subroutine init_equilibrium(f1,r,ux,uy)
  use declare
  real :: cidotu
  integer :: x,y,i
  real :: f1(1:nx,1:ny,0:ndir)
  real :: r(1:nx,1:ny),ux(1:nx,1:ny),uy(1:nx,1:ny)
  DO x= 1,nx
    DO y =1,ny
      DO i = 0,ndir
        cidotu = dirx(i)*ux(x,y)+diry(i)*uy(x,y)
        f1(x,y,i) = wi(i)*r(x,y)*(1+3*cidotu+4.5*cidotu*cidotu-1.5*(ux(x,y)*ux(x,y)+uy(x,y)*uy(x,y)))
      END DO
    END DO
  END DO
  !print*,f1(1,1,0)
END subroutine init_equilibrium


subroutine stream(r,f_src,f_dst)
  use declare
  real :: f_src(1:nx,1:ny,0:ndir),f_dst(1:nx,1:ny,0:ndir)
  integer :: x,y,xmd,ymd,i
  real :: r(1:nx,1:ny)
  DO x = 1,nx
    DO y = 2,ny-1
      DO i = 0,ndir
        xmd = MOD(NX+x-dirx(i)-1,NX)+1
        ymd = MOD(NY+y-diry(i)-1,NY)+1
        f_dst(x,y,i) = f_src(xmd,ymd,i)
        !print *, x,y,i,xmd,ymd,f_src(xmd,ymd,i),f_dst(xmd,ymd,i)
      ENDDO
    ENDDO
    y = 1
    DO i = 0,ndir
      IF((i.EQ.0).OR.(i.EQ.1).OR.(i.EQ.3).OR.(i.EQ.4).OR.(i.EQ.7).OR.(i.EQ.8)) THEN
        xmd = MOD(NX+x-dirx(i)-1,NX)+1
        ymd = MOD(NY+y-diry(i)-1,NY)+1
        f_dst(x,y,i) = f_src(xmd,ymd,i)
      ELSEIF((i.EQ.2).OR.(i.EQ.5).OR.(i.EQ.6)) THEN
        cidotu = dirx(pn(i))*uLowWallx+diry(pn(i))*uLowWally
        f_dst(x,y,i) = f_src(x,y,pn(i))-6*wi(pn(i))*r(x,y)*cidotu
      ENDIF
    ENDDO
    y = NY
    DO i = 0,ndir
      IF((i.EQ.0).OR.(i.EQ.1).OR.(i.EQ.2).OR.(i.EQ.3).OR.(i.EQ.5).OR.(i.EQ.6)) THEN
        xmd = MOD(NX+x-dirx(i)-1,NX)+1
        ymd = MOD(NY+y-diry(i)-1,NY)+1
        f_dst(x,y,i) = f_src(xmd,ymd,i)
      ELSEIF((i.EQ.4).OR.(i.EQ.7).OR.(i.EQ.8)) THEN
        cidotu = dirx(pn(i))*uUpWallx+diry(pn(i))*uUpWally
        f_dst(x,y,i) = f_src(x,y,pn(i))
      ENDIF
    END DO
  END DO
END subroutine stream


subroutine compute_rho_u(f2,r,ux,uy)
  use declare
  real :: f2(1:nx,1:ny,0:ndir)
  real :: r(1:nx,1:ny),ux(1:nx,1:ny),uy(1:nx,1:ny)
  integer :: x,y,i
  DO x = 1,nx
    DO y = 1,ny
      r(x,y) = 0
      ux(x,y) = 0
      uy(x,y) = 0
      DO i = 0,ndir
        r(x,y)  = r(x,y)  + f2(x,y,i)
        ux(x,y) = ux(x,y) + f2(x,y,i)*dirx(i)
        uy(x,y) = uy(x,y) + f2(x,y,i)*diry(i)
      END DO
      ux(x,y) = ux(x,y)/r(x,y)
      uy(x,y) = uy(x,y)/r(x,y)
    END DO
  END DO
END subroutine compute_rho_u


subroutine collide(f,r,ux,uy,f_out)
  use declare
  real :: f(1:nx,1:ny,0:ndir),f_out(1:nx,1:ny,0:ndir),feq(1:nx,1:ny,0:ndir)
  real :: r(1:nx,1:ny),ux(1:nx,1:ny),uy(1:nx,1:ny)
  integer :: x,y,i
  real :: cidotu
  DO x = 1,nx
    DO y =1,ny
      DO i = 0,ndir
        cidotu = dirx(i)*ux(x,y)+diry(i)*uy(x,y)
        feq(x,y,i) = wi(i)*r(x,y)*(1+3*cidotu+4.5*cidotu*cidotu-1.5*(ux(x,y)*ux(x,y)+uy(x,y)*uy(x,y)))
        f_out(x,y,i) = (omtauniv*f(x,y,i))+(tauinv*feq(x,y,i))
      END DO
    END DO
  END DO
END subroutine collide
