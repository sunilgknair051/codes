module declare
    implicit none
    real :: w0,ws,wd
    parameter(w0=4./9., ws=1./9., wd= 1./36.)
    integer :: nx,ny,ndir
    parameter(ny=200,nx=150,ndir=8)
    integer :: dirx(0:8)=(/0,1,0,-1,0,1,-1,-1,1/)
    integer :: diry(0:8)=(/0,0,1,0,-1,1,1,-1,-1/)
    integer, parameter :: no_slip(0:8)=(/0,3,4,1,2,7,8,5,6/)
    real    :: wi(0:8)=(/w0,ws,ws,ws,ws,wd,wd,wd,wd/)
    real,parameter :: nu = 1./6.
    real,parameter :: rho0 = 1.
    real,parameter :: tauinv = 2./(6.*nu + 1.)
    real,parameter :: omtauinv = 1. - tauinv
    real,parameter :: scale = 1.
    real,parameter :: x_sphere = nx/2., y_sphere = 1., rad_sphere = 0., PI = 4.D0*DATAN(1.D0)
    real,parameter :: delfac = 1./6., delta = ny*delfac, omega = 2.*nu/delta**2.
    !integerpy, parameter :: NSTEPS = int(3./omega)
    integer, parameter :: NSTEPS = 500
    real,parameter :: Zq_SI = 8.8E6, nu_SI = 1E-6, rho_Liq_SI = 1000 , f0_SI = 5E6, delta_SI = 250E-9
    real,parameter :: uxWall = 1.
    integer:: beginning, rate, end
end module declare


program repeat
!!$OMP PARALLEL
use declare
use omp_lib
implicit none
integer :: x,y,i
complex :: dfc
complex :: dr(1:nx,1:ny),ux(1:nx,1:ny),uy(1:nx,1:ny)
complex :: FxLiqOnSphe, FxLiqOnWall
complex :: h1(1:nx,1:ny,0:ndir),h2(1:nx,1:ny,0:ndir)


DO x = 1,nx
  DO y = 1,ny
    dr(x,y) = complex(0,0)
    ux(x,y) = complex(0,0)
    uy(x,y) = complex(0,0)
  ENDDO
ENDDO


open(unit = 20,file='complex_freq2.asc',action='write')

call init_equilibrium(h1,dr,ux,uy)
write (20,*) NSTEPS
call system_clock(beginning,rate)
!$OMP PARALLEL do
do i = 1,NSTEPS
  call stream(h1,h2,FxLiqOnSphe, FxLiqOnWall)
  call compute_dr_u(h2,dr,ux,uy)
  call collide(h2,dr,ux,uy,h1)
  call calc_dfc(FxLiqOnSphe, FxLiqOnWall,dfc)
  write (20,*) real(dfc),aimag(dfc)
  print *,real(dfc), aimag(dfc)
  !print *, "Time = ",ETIME
end do
!$OMP END PARALLEL DO
close(unit = 20)

open(unit = 30,file='uxuy2.asc',action='write')
write(30,*)nx,ny
DO x=1,nx
  DO y=1,ny
    write(30,*)x,y,ux(x,y),uy(x,y)
  ENDDO
ENDDO
close(unit = 30)

!!$OMP END PARALLEL
call system_clock(end)
print*,"Elapsed time is in seconds --->",real(end-beginning)/real(rate)
end program repeat

subroutine init_equilibrium(h1,dr,ux,uy)
  use declare
  implicit none
  complex :: cidotu
  integer :: x,y,i
  complex :: h1(1:nx,1:ny,0:ndir)
  complex :: dr(1:nx,1:ny),ux(1:nx,1:ny),uy(1:nx,1:ny)
  !$OMP PARALLEL do
  DO x= 1,nx
    DO y =1,ny
      DO i = 0,ndir
        cidotu = dirx(i)*ux(x,y)+diry(i)*uy(x,y)
        h1(x,y,i) = wi(i)*dr(x,y)*(1+3*cidotu)
      END DO
    END DO
  END DO
  !print*,f1(1,1,0)
  !$OMP END PARALLEL DO
END subroutine init_equilibrium


subroutine stream(h_src,h_dst,FxLiqOnSphe,FxLiqOnWall)
  use declare
  implicit none
  complex :: h_src(1:nx,1:ny,0:ndir),h_dst(1:nx,1:ny,0:ndir)
  integer :: x,y,xmd,ymd,i,ibar
  complex :: FxLiqOnWall, FxLiqOnSphe
  FxLiqOnWall = 0
  FxLiqOnSphe = 0
  !!$OMP PARALLEL do
  DO x = 1,nx
    DO y = 1,ny
      DO i = 0,ndir
        IF ((y.le.ny).and.(y.ge.1).and.((x-x_sphere)**2 +(y-y_sphere)**2.gt.rad_sphere**2)) THEN
          ibar = no_slip(i)
          xmd = MOD(NX+x-dirx(i)-1,NX)+1
          ymd = y-diry(i)
          IF((ymd.le.ny).and.(ymd.ge.1).and.((xmd-x_sphere)**2 +(ymd-y_sphere)**2.ge.rad_sphere**2)) THEN
            h_dst(x,y,i) = h_src(xmd,ymd,i)
          ELSEIF(ymd.gt.ny) THEN
            h_dst(x,y,i) = h_src(x,y,ibar)
          ELSEIF(ymd.lt.1) THEN
            h_dst(x,y,i) = h_src(x,y,ibar) - 6*wi(ibar)*dirx(ibar)*uxWall
            FxLiqOnWall = FxLiqOnWall + (h_dst(x,y,i)+h_src(x,y,ibar))*dirx(i)
          ELSEIF((xmd-x_sphere)**2 +(ymd-y_sphere)**2.lt.rad_sphere**2) THEN
            h_dst(x,y,i) = h_src(x,y,ibar) - 6*wi(ibar)*dirx(ibar)*uxWall
            FxLiqOnSphe = FxLiqOnSphe + (h_dst(x,y,i)+h_src(x,y,ibar))*dirx(i)
          ENDIF
        ENDIF
        !print *, x,y,i,xmd,ymd,f_src(xmd,ymd,i),f_dst(xmd,ymd,i)
      END DO
    END DO
  END DO
  !!$OMP END PARALLEL do
END subroutine stream


subroutine compute_dr_u(h_in,dr,ux,uy) !h_in...h2
  use declare
  implicit none
  complex :: h_in(1:nx,1:ny,0:ndir)
  complex :: dr(1:nx,1:ny),ux(1:nx,1:ny),uy(1:nx,1:ny)
  integer :: x,y,i
  !$OMP PARALLEL do
  DO x = 1,nx
    DO y = 1,ny
      dr(x,y) = complex(0,0)
      ux(x,y) = complex(0,0)
      uy(x,y) = complex(0,0)
      IF ((y.le.ny).and.(y.ge.1).and.((x-x_sphere)**2 +(y-y_sphere)**2.ge.rad_sphere**2)) THEN
        DO i = 0,ndir
          dr(x,y) = dr(x,y) + h_in(x,y,i)
          ux(x,y) = ux(x,y) + h_in(x,y,i)*dirx(i)
          uy(x,y) = uy(x,y) + h_in(x,y,i)*diry(i)
        ENDDO
      ENDIF
    ENDDO
  ENDDO
  !$OMP END PARALLEL do
END subroutine compute_dr_u


subroutine collide(h_in,dr,ux,uy,h_out)
  use declare
  implicit none
  complex,intent(in) :: h_in(1:nx,1:ny,0:ndir)
  complex,intent(out) :: h_out(1:nx,1:ny,0:ndir)
  complex :: heq
  complex,intent(in) :: dr(1:nx,1:ny),ux(1:nx,1:ny),uy(1:nx,1:ny)
  integer :: x,y,i
  real, parameter:: vomega = omega**2/2
  complex :: cidotu
  !$OMP PARALLEL do
  DO x = 1,nx
    DO y =1,ny
      IF ((y.le.ny).and.(y.ge.1).and.((x-x_sphere)**2 +(y-y_sphere)**2.ge.rad_sphere**2)) THEN
        DO i = 0,ndir
          cidotu = dirx(i)*ux(x,y) + diry(i)*uy(x,y)
          heq = wi(i) * (dr(x,y) + 3*cidotu)
          h_out(x,y,i) = omtauinv*h_in(x,y,i) + tauinv*heq +complex(0,1)*omega*h_in(x,y,i) - vomega*heq
        END DO
      ENDIF
    ENDDO
  ENDDO
  !$OMP END PARALLEL do
END subroutine collide

subroutine calc_dfc(FxLiqOnSphe, FxLiqOnWall,dfc)
  use declare
  implicit none
  complex :: TanStress, FxLiqOnSphe, FxLiqOnWall, ZL, Zq, dfc
  TanStress = (FxLiqOnWall+FxLiqOnSphe)/nx
  ZL = TanStress / uxWall
  Zq = (rho0/rho_Liq_SI)*(nu/nu_SI)*(delta_SI/delta)*Zq_SI
  dfc = f0_SI*complex(0,1)/PI *ZL/Zq
  dfc = -CONJG(dfc)
end subroutine calc_dfc
