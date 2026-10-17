
module declare
    implicit none
    real*8 :: w0,ws,wd
    parameter(w0=1./3., ws=1./18., wd= 1./36.)
    integer :: nx,ny,nz,ndir
    integer,parameter :: ncycles = 10
    parameter(nx=1,ny=250,nz=1,ndir=18)
    integer :: dirx(0:18)=(/0,1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0/)
    integer :: diry(0:18)=(/0,0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1/)
    integer :: dirz(0:18)=(/0,0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1/)
    integer, parameter :: no_slip(0:18)=(/0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17/)
    real*8   :: wi(0:18)=(/w0,ws,ws,ws,ws,ws,ws,wd,wd,wd,wd,wd,wd,wd,wd,wd,wd,wd,wd/)
    real*8,parameter :: nu = 1./6.
    real*8,parameter :: rho0 = 1.
    real*8,parameter :: tauinv = 2./(6.*nu + 1.)
    real*8,parameter :: omtauinv = 1. - tauinv
    real*8,parameter :: scale = 1.
    real*8,parameter :: x_sphere = nx/2., y_sphere = 0., z_sphere = nz/2., rad_sphere = 0., PI = 4.D0*DATAN(1.D0)
    real*8,parameter :: delfac = 1./3., delta = ny*delfac, omega = 2.*nu/delta**2.
    integer, parameter :: NSTEPS = int(ncycles*2*PI/omega)
    !integer, parameter :: Nsteps = 2
    ! integer, parameter :: NSTEPS = 50000
    real*8,parameter :: Zq_SI = 8.8E6, nu_SI = 1E-6, rho_Liq_SI = 1000 , f0_SI = 5E6
    real*8,parameter :: delta_SI = sqrt(2*nu_SI/(2*PI*f0_SI))
    real*8,parameter :: uxWall = 1.
    integer:: beginning, rate, end
end module declare


program repeat
!!$OMP PARALLEL
use declare
use omp_lib
implicit none
integer :: x,y,z,i
double complex :: dfc
double complex :: dr(1:nx,1:ny,1:nz),ux(1:nx,1:ny,1:nz),uy(1:nx,1:ny,1:nz),uz(1:nx,1:ny,1:nz)
double complex :: FxLiqOnSphe, FxLiqOnWall
double complex :: h1(1:nx,1:ny,1:nz,0:ndir),h2(1:nx,1:ny,1:nz,0:ndir)

DO x = 1,nx
  DO y = 1,ny
    DO z = 1,nz
     dr(x,y,z) = complex(0,0)
     ux(x,y,z) = complex(0,0)
     uy(x,y,z) = complex(0,0)
     uz(x,y,z) = complex(0,0)
    ENDDO
  ENDDO
ENDDO

open(unit = 20,file='3d_dp_complex_freq.asc',action='write')

call init_equilibrium(h1,dr,ux,uy,uz)
write (20,*) NSTEPS
call system_clock(beginning,rate)
!!$OMP PARALLEL do num_threads(8)
do i = 1,NSTEPS
  call stream(h1,h2,FxLiqOnSphe, FxLiqOnWall)
  call compute_dr_u(h2,dr,ux,uy,uz)
  call collide(h2,dr,ux,uy,uz,h1)
  call calc_dfc(FxLiqOnSphe, FxLiqOnWall,dfc)
  write (20,*) real(dfc),aimag(dfc)
  if (mod(1,10).eq.0) then
    print *,1,NSTEPS,real(dfc), aimag(dfc)
  endif
  !print *, "Time = ",ETIME
end do
!!$OMP END PARALLEL do
close(unit = 20)

open(unit = 30,file='3d_dp_uxuy.asc',action='write')
write(30,*)nx,ny,nz
DO x=1,nx
  DO y=1,ny
    DO z=1,nz
      write(30,*)x,y,z,ux(x,y,z),uy(x,y,z),uz(x,y,z)
    ENDDO
  ENDDO
ENDDO
close(unit = 30)

!!$OMP END PARALLEL
call system_clock(end)
print*,"Elapsed time is in seconds --->",real(end-beginning)/real(rate)
end program repeat

subroutine init_equilibrium(h1,dr,ux,uy,uz)
  use declare
  implicit none
  double complex :: cidotu
  integer :: x,y,z,i
  double complex :: h1(1:nx,1:ny,1:nz,0:ndir)
  double complex :: dr(1:nx,1:ny,1:nz),ux(1:nx,1:ny,1:nz),uy(1:nx,1:ny,1:nz),uz(1:nx,1:ny,1:nz)
!!$OMP PARALLEL do num_threads(8)
  DO x= 1,nx
    DO y =1,ny
      DO z =1,nz
        DO i = 0,ndir
         cidotu = dirx(i)*ux(x,y,z)+diry(i)*uy(x,y,z)+dirz(i)*uz(x,y,z)
         h1(x,y,z,i) = wi(i)*dr(x,y,z)*(1+3*cidotu)
        END DO
      ENDDO
    END DO
  END DO
  !print*,f1(1,1,0)
  !!$OMP END PARALLEL DO
END subroutine init_equilibrium


subroutine stream(h_src,h_dst,FxLiqOnSphe,FxLiqOnWall)
  use declare
  implicit none
  double complex :: h_src(1:nx,1:ny,1:nz,0:ndir),h_dst(1:nx,1:ny,1:nz,0:ndir)
  integer :: x,y,z,xmd,ymd,zmd,i,ibar
  double complex :: FxLiqOnWall, FxLiqOnSphe
  FxLiqOnWall = 0
  FxLiqOnSphe = 0
  !!$OMP PARALLEL do num_threads(8)
  DO x = 1,nx
    DO y = 1,ny
      DO z = 1,nz
        DO i = 0,ndir
         IF ((y.le.ny).and.(y.ge.1).and.((x-x_sphere)**2 +(y-y_sphere)**2+(z-z_sphere)**2.gt.rad_sphere**2)) THEN
           ibar = no_slip(i)
           xmd = MOD(NX+x-dirx(i)-1,NX)+1
           ymd = y-diry(i)
           zmd = MOD(NZ+z-dirz(i)-1,NZ)+1
            IF((ymd.le.ny).and.(ymd.ge.1).and.((xmd-x_sphere)**2 +(ymd-y_sphere)**2+(zmd-z_sphere)**2.ge.rad_sphere**2)) THEN
              h_dst(x,y,z,i) = h_src(xmd,ymd,zmd,i)
            ELSEIF(ymd.gt.ny) THEN
              h_dst(x,y,z,i) = h_src(x,y,z,ibar)
            ELSEIF(ymd.lt.1) THEN
              h_dst(x,y,z,i) = h_src(x,y,z,ibar) - 6*wi(ibar)*dirx(ibar)*uxWall
              FxLiqOnWall = FxLiqOnWall + (h_dst(x,y,z,i)+h_src(x,y,z,ibar))*dirx(i)
            ELSEIF((xmd-x_sphere)**2 +(ymd-y_sphere)**2+(zmd-z_sphere)**2.lt.rad_sphere**2) THEN
              h_dst(x,y,z,i) = h_src(x,y,z,ibar) - 6*wi(ibar)*dirx(ibar)*uxWall
              FxLiqOnSphe = FxLiqOnSphe + (h_dst(x,y,z,i)+h_src(x,y,z,ibar))*dirx(i)
            ENDIF
          ENDIF
              !print *, x,y,i,xmd,ymd,f_src(xmd,ymd,i),f_dst(xmd,ymd,i)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !!$OMP END PARALLEL do
END subroutine stream


subroutine compute_dr_u(h_in,dr,ux,uy,uz) !h_in...h2
  use declare
  implicit none
  double complex :: h_in(1:nx,1:ny,1:nz,0:ndir)
  double complex :: dr(1:nx,1:ny,1:nz),ux(1:nx,1:ny,1:nz),uy(1:nx,1:ny,1:nz),uz(1:nx,1:ny,1:nz)
  integer :: x,y,z,i
  !!$OMP PARALLEL do  num_threads(8)
  DO x = 1,nx
    DO y = 1,ny
      DO z = 1,nz
        dr(x,y,z) = complex(0,0)
        ux(x,y,z) = complex(0,0)
        uy(x,y,z) = complex(0,0)
        uz(x,y,z) = complex(0,0)
        IF ((y.le.ny).and.(y.ge.1).and.((x-x_sphere)**2 +(y-y_sphere)**2+(z-z_sphere)**2.ge.rad_sphere**2)) THEN
          DO i = 0,ndir
           dr(x,y,z) = dr(x,y,z) + h_in(x,y,z,i)
           ux(x,y,z) = ux(x,y,z) + h_in(x,y,z,i)*dirx(i)
           uy(x,y,z) = uy(x,y,z) + h_in(x,y,z,i)*diry(i)
           uz(x,y,z) = uz(x,y,z) + h_in(x,y,z,i)*dirz(i)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  !!$OMP END PARALLEL do
END subroutine compute_dr_u


subroutine collide(h_in,dr,ux,uy,uz,h_out)
  use declare
  implicit none
  double complex,intent(in)  :: h_in( 1:nx,1:ny,1:nz,0:ndir)
  double complex,intent(out) :: h_out(1:nx,1:ny,1:nz,0:ndir)
  double complex :: heq
  double complex,intent(in) :: dr(1:nx,1:ny,1:nz),ux(1:nx,1:ny,1:nz),uy(1:nx,1:ny,1:nz),uz(1:nx,1:ny,1:nz)
  integer :: x,y,z,i
  real, parameter:: vomega = omega**2/2
  double complex :: cidotu
  !!$OMP PARALLEL do num_threads(8)
  DO x = 1,nx
    DO y = 1,ny
      DO z = 1,nz
        IF ((y.le.ny).and.(y.ge.1).and.((x-x_sphere)**2 +(y-y_sphere)**2+(z-z_sphere)**2.ge.rad_sphere**2)) THEN
          DO i = 0,ndir
            cidotu = dirx(i)*ux(x,y,z) + diry(i)*uy(x,y,z)+ dirz(i)*uz(x,y,z)
            heq = wi(i) * (dr(x,y,z) + 3*cidotu)
            h_out(x,y,z,i) = omtauinv*h_in(x,y,z,i) + tauinv*heq +complex(0,1)*omega*h_in(x,y,z,i) - vomega*heq
          END DO
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  !!$OMP END PARALLEL do
END subroutine collide

subroutine calc_dfc(FxLiqOnSphe, FxLiqOnWall,dfc)
  use declare
  implicit none
  double complex :: TanStress, FxLiqOnSphe, FxLiqOnWall, ZL, Zq, dfc
  TanStress = (FxLiqOnWall+FxLiqOnSphe)/nx
  ZL = TanStress / uxWall
  Zq = (rho0/rho_Liq_SI)*(nu/nu_SI)*(delta_SI/delta)*Zq_SI
  dfc = f0_SI*complex(0,1)/PI *ZL/Zq
  dfc = -CONJG(dfc)
end subroutine calc_dfc
