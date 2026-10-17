module declare
  implicit none
  integer, parameter :: ndir = 18, nx_mx = 100, nz_mx = 100, ny_mx = 100
  integer, parameter :: dirx(0:18) = (/0,+1,-1, 0, 0, 0, 0,+1,-1,+1,-1, 0, 0,+1,-1,+1,-1, 0, 0/)
  integer, parameter :: dirz(0:18) = (/0, 0, 0, 0, 0,+1,-1, 0, 0,+1,-1,+1,-1, 0, 0,-1,+1,-1,+1/)
  integer, parameter :: diry(0:18) = (/0, 0, 0,+1,-1, 0, 0,+1,-1, 0, 0,+1,-1,-1,+1, 0, 0,+1,-1/)
  integer, parameter :: no_slip(0:18) = (/0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17/)
  real*8,parameter :: pi = 4.D0*DATAN(1.D0), nperiods = 1, nu = 1./6.,y_sphere = +1.
  real*8,parameter :: wi(0:18) = (/ 1./3.,1./18,1./18.,1./18.,1./18.,1./18.,1./18., &
                           1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36./)
  real*8,parameter :: rho0 = 1., tauinv = 2./(6.*nu + 1.), omtauinv = 1. - tauinv
  real*8,parameter :: nu_SI = 1e-6, rho_Liq_SI = 1000., f0_SI = 5e6, delta_f0_SI = sqrt(2 * nu_SI/(2*pi*f0_SI));
  real*8,parameter :: rho_sphere_SI = 1.e3
  real*8           :: rad_sphere_SI, cell_height0_SI, cell_height_SI, gridreso0_SI, gridreso_SI

  integer :: nx,ny,nz,nsize,nreso,nsteps
  real*8 :: x_sphere,z_sphere,rad_sphere,delta,omega,omega2by2
  complex*16 :: dr(1:nx_mx,1:ny_mx,1:nz_mx), ux(1:nx_mx,1:ny_mx,1:nz_mx), uy(1:nx_mx,1:ny_mx,1:nz_mx), uz(1:nx_mx,1:ny_mx,1:nz_mx)
  complex*16 :: h1(1:nx_mx,1:ny_mx,1:nz_mx,0:ndir), h2(1:nx_mx,1:ny_mx,1:nz_mx,0:ndir)
  complex*16 :: FxLiqOnSphe,FOS,ratio
end module declare

program FT_LBM
  use declare; implicit none
  integer :: it,ireso,isize

  nreso = 1; nsize = 1;
!  rad_sphere_SI = 100e-9; cell_height0_SI = 1200e-9; gridreso0_SI = 15e-9; -> 0.8, 0.6  5.1.2019
  rad_sphere_SI = 100e-9; cell_height0_SI = 1200e-9; gridreso0_SI = 15e-9; !-> 0.8, 0.6

  do ireso = 1,nreso
    gridreso_SI = gridreso0_SI! * (1. + 1.*(ireso-1)/(nreso-1));
    delta = delta_f0_SI/gridreso_SI
    rad_sphere = rad_sphere_SI/gridreso_SI;
    do isize  = 1,nsize
      cell_height_SI = cell_height0_SI! * (1. + 1.*(isize-1)/(nsize-1));
      ny = int(cell_height_SI/gridreso_SI); nx = ny; nz = ny
      x_sphere = nx/2.; z_sphere = nz/2.! y_sphere = ny/2.;
      omega = 2.*nu/delta**2.; omega2by2 = omega**2/2;
      nsteps = int(3./omega);
      print *,omega,nsteps,ny
      call init_u0;
      do it = 1,nsteps
        call LB_Stream; call compute_dr_u; call collide; call calc_dfc;
        if(mod(it,20) == 1) print *,it,nsteps,real(ratio),aimag(ratio)
      enddo
    enddo
  enddo
end program FT_LBM

subroutine init_u0
  use declare; implicit none
  integer :: x,y,z,i
  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        do i = 0,ndir
          h1(x,y,z,i) = 0
        end do
      end do
    end do
  end do
end subroutine init_u0

subroutine LB_Stream
  use declare; implicit none
  integer :: x,y,z,xmd,ymd,zmd,i,ibar
  FxLiqOnSphe = 0;
  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        do i = 0,ndir
          if ((x-x_sphere)**2 +(y-y_sphere)**2 +(z-z_sphere)**2.gt.rad_sphere**2) then
            ibar = no_slip(i)
            xmd = mod(nx+x-dirx(i)-1,nx)+1
            ymd = mod(ny+y-diry(i)-1,ny)+1
            zmd = mod(nz+z-dirz(i)-1,nz)+1
            if((xmd-x_sphere)**2 +(ymd-y_sphere)**2 +(zmd-z_sphere)**2.ge.rad_sphere**2) then
              h2(x,y,z,i) = h1(xmd,ymd,zmd,i)
            elseif((xmd-x_sphere)**2 + (ymd-y_sphere)**2 + (zmd-z_sphere)**2.lt.rad_sphere**2) then
              h2(x,y,z,i) = h1(x,y,z,ibar) - 6*wi(ibar)*dirx(ibar)
              FxLiqOnSphe = FxLiqOnSphe + (h2(x,y,z,i) + h1(x,y,z,ibar))*dirx(i)
            endif
          endif
        end do
      end do
    end do
  end do
end subroutine LB_Stream

subroutine compute_dr_u
  use declare; implicit none
  integer :: x,y,z,i
  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        if ((y.le.ny).and.(y.ge.1).and.((x-x_sphere)**2 +(y-y_sphere)**2 +(z-z_sphere)**2.ge.rad_sphere**2)) then
          dr(x,y,z) = complex(0,0)
          ux(x,y,z) = complex(0,0)
          uy(x,y,z) = complex(0,0)
          uz(x,y,z) = complex(0,0)
          do i = 0,ndir
            dr(x,y,z) = dr(x,y,z) + h2(x,y,z,i)
            ux(x,y,z) = ux(x,y,z) + h2(x,y,z,i)*dirx(i)
            uy(x,y,z) = uy(x,y,z) + h2(x,y,z,i)*diry(i)
            uz(x,y,z) = uz(x,y,z) + h2(x,y,z,i)*dirz(i)
          enddo
        endif
      enddo
    enddo
  enddo
end subroutine compute_dr_u

subroutine collide
  use declare; implicit none
  complex*16 :: heq,cidotu
  integer :: x,y,z,i
  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        if ((y.le.ny).and.(y.ge.1).and.((x-x_sphere)**2 +(y-y_sphere)**2 +(z-z_sphere)**2.ge.rad_sphere**2)) then
          do i = 0,ndir
            cidotu = dirx(i)*ux(x,y,z) + diry(i)*uy(x,y,z) + dirz(i)*uz(x,y,z)
            heq = wi(i) * (dr(x,y,z) + 3*cidotu)
!            h1(x,y,z,i) = omtauinv * h2(x,y,z,i) + tauinv*heq +complex(0,1)*omega*h2(x,y,z,i) - omega**2/2*heq
            h1(x,y,z,i) = tauinv*heq + complex(0,1)*omega*h2(x,y,z,i) - omega2by2 *heq  ! omtauinv is zero
          end do
        endif
      enddo
    enddo
  enddo
end subroutine collide

subroutine calc_dfc
  use declare; implicit none
  real*8 :: FOS_analy_Re,FOS_analy_Im, eta
  eta = nu * rho0
  FOS = FxLiqOnSphe;
  FOS_analy_Re = 6. * pi * eta * rad_sphere *(1. + rad_sphere/delta)
  FOS_analy_Im = 4.*pi/3.*rho0*rad_sphere**3 * (1./2. + 9./4.*delta/rad_sphere)
  ratio = (FOS_analy_Re + complex(0,1)*omega * FOS_analy_Im) / FOS

end subroutine calc_dfc
