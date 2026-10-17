module declare
  implicit none
  integer, parameter :: ndir = 18, nx_mx = 70, nz_mx = 70, ny_mx = 50, novt_mx = 5
  integer, parameter :: dirx(0:18) = (/0,+1,-1, 0, 0, 0, 0,+1,-1,+1,-1, 0, 0,+1,-1,+1,-1, 0, 0/)
  integer, parameter :: dirz(0:18) = (/0, 0, 0, 0, 0,+1,-1, 0, 0,+1,-1,+1,-1, 0, 0,-1,+1,-1,+1/)
  integer, parameter :: diry(0:18) = (/0, 0, 0,+1,-1, 0, 0,+1,-1, 0, 0,+1,-1,-1,+1, 0, 0,+1,-1/)
  integer, parameter :: no_slip(0:18) = (/0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17/)
  real*8,parameter :: pi = 4.D0*DATAN(1.D0), y_sphere = +5., nperiods = 1, nu = 1./6.
  real*8,parameter :: wi(0:18) = (/ 1./3.,1./18,1./18.,1./18.,1./18.,1./18.,1./18., &
                           1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36./)
  real*8,parameter :: rho0 = 1., tauinv = 2./(6.*nu + 1.), omtauinv = 1. - tauinv, uxWall = 1., delfac = 1.
  real*8,parameter :: Zq_SI = 8.8e6, nu_SI = 1e-6, rho_Liq_SI = 1000., f0_SI = 5e6, delta_f0_SI = sqrt(2 * nu_SI/(2*pi*f0_SI));
  real*8,parameter :: rho_sphere_SI = 1.e3
  real*8           :: rad_sphere_SI, cell_height_SI, gridreso_SI, delta_f0

  integer :: beginning, rate, ending,nx,ny,nz,nsize,novt,nsteps,ovtorder
  real*8 :: x_sphere,z_sphere,rad_sphere,Sizefac,delta,omega,omega2by2,f_SI
  complex*16 :: dr(1:nx_mx,1:ny_mx,1:nz_mx), ux(1:nx_mx,1:ny_mx,1:nz_mx), uy(1:nx_mx,1:ny_mx,1:nz_mx), uz(1:nx_mx,1:ny_mx,1:nz_mx)
  complex*16 :: h1(1:nx_mx,1:ny_mx,1:nz_mx,0:ndir), h2(1:nx_mx,1:ny_mx,1:nz_mx,0:ndir)
  complex*16 :: FxLiqOnSphe, FxLiqOnWall, dfc,TanStress,ZL,Zq,uxTop,uyTop,Zliq
  complex*16 :: dfc_ref(1:novt_mx),dfc_fin(1:novt_mx)
end module declare

program FT_LBM
  use declare; implicit none
  integer :: x,y,z,it,isize,iovt
  open(unit = 20,file='complex_freq.asc1',action='write',status='unknown'); close(unit = 20,status='delete');
  open(unit = 21,file='fin_freq.asc1',    action='write',status='unknown'); close(unit = 21,status='delete');


  nsize = 20; novt = 2; nsteps = 10001
  rad_sphere_SI = 5e-9; cell_height_SI = 15e-9; gridreso_SI = 0.1e-9;
  ny = int(cell_height_SI/gridreso_SI)
  delta_f0 = delta_f0_SI * float(ny)/cell_height_SI

  do iovt = 1,novt
    nx = 1; nz = 1; rad_sphere = 0;
    ovtorder = 2 * iovt - 1; f_SI = ovtorder* f0_SI;
    delta = delta_f0 / sqrt(float(ovtorder));
    omega = 2.*nu/delta**2.; omega2by2 = omega**2/2; Zliq = -conjg(sqrt(rho0*complex(0,1)*omega*rho0*nu))
    call init_u0;
    do it = 1,nsteps
      call LB_Stream; call compute_dr_u; call collide; call calc_dfc;
!      if ((it.ge.1500).and.(mod(it,100)==0)) then
      if ((it.ge.500)) then
        open(unit = 20,file='complex_freq.asc1',action='write',status='unknown', position='append');
        write (20,*) isize,nsize,sizefac,iovt,it,real(dfc),aimag(dfc);close(unit = 20);
      endif
      if(mod(it,500) == 1) print *,it,real(dfc),aimag(dfc)
    enddo
    open(unit = 21,file='fin_freq.asc1',    action='write',status='unknown', position='append');
    write (21,*) isize,Sizefac,iovt,real(dfc),aimag(dfc); close(unit = 21)
    dfc_ref(iovt) = dfc
  enddo

  rad_sphere = rad_sphere_SI * delta_f0 / delta_f0_SI;
  print *,'rad_sphere , delta, ny',rad_sphere,delta,ny
  do isize = 1,nsize
    call system_clock(beginning,rate)
    if (nsize.gt.1) Sizefac = 1.2 + (isize-1.) / (nsize-1.)*2.
    if (nsize.eq.1) Sizefac = 1.2
    rad_sphere = 10.; nx = int(2 * rad_sphere * Sizefac); nz = nx; x_sphere = nx/2.; z_sphere = nz/2.
    do iovt = 1,novt
      ovtorder = 2 * iovt - 1;
      delta = delta_f0 / sqrt(float(ovtorder));
      omega = 2.*nu/delta**2.; omega2by2 = omega**2/2; Zliq = -conjg(sqrt(rho0*complex(0,1)*omega*rho0*nu))
      call init_u0;
      do it = 1,nsteps
        call LB_Stream; call compute_dr_u; call collide; call calc_dfc;
        if (it.ge.500) then
          open(unit = 20,file='complex_freq.asc1',action='write',status='unknown', position='append');
          write (20,*) it,real(dfc),aimag(dfc);close(unit = 20);
        endif
        if(mod(it,500) == 1) print *,it,nsteps,real(dfc),aimag(dfc)
      enddo
      dfc_fin(iovt) = dfc - dfc_ref(iovt)
      print *,'iovt',iovt,dfc,dfc_ref(iovt),dfc_fin(iovt)
      print *,'fin',isize,nsize,ovtorder,real(dfc_fin/ovtorder),aimag(dfc_fin/ovtorder)
      call system_clock(ending); print *,"Elapsed time is in seconds --->",real(ending - beginning)/real(rate)
    enddo
    open(unit = 21,file='fin_freq.asc1',    action='write',status='unknown', position='append');
    write (21,*) isize,pi/4./Sizefac**2,real(dfc_fin)/ovtorder,aimag(dfc_fin)/ovtorder; close(unit = 21)
  enddo
  open(unit = 30,file='uxuy_3D.asc1',action='write')
  write(30,"(I5,I5,I5)")nx,ny,nz
  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        write(30,"(I5,I5,I5,F12.8,F12.8,F12.8,F12.8,F12.8,F12.8)")x,y,z,real(ux(x,y,z)),aimag(ux(x,y,z)),&
                                       real(uy(x,y,z)),aimag(uy(x,y,z)),real(uz(x,y,z)),aimag(uz(x,y,z))
      enddo
    enddo
  enddo
  close(unit = 30)
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
  FxLiqOnWall = 0; FxLiqOnSphe = 0;
  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        do i = 0,ndir
          if ((y.le.ny).and.(y.ge.1).and.((x-x_sphere)**2 +(y-y_sphere)**2 +(z-z_sphere)**2.gt.rad_sphere**2)) then
            ibar = no_slip(i)
            xmd = mod(nx+x-dirx(i)-1,nx)+1
            ymd = y-diry(i)
            zmd = mod(nz+z-dirz(i)-1,nz)+1
            if((ymd.le.ny).and.(ymd.ge.1).and.((xmd-x_sphere)**2 +(ymd-y_sphere)**2 +(zmd-z_sphere)**2.ge.rad_sphere**2)) then
              h2(x,y,z,i) = h1(xmd,ymd,zmd,i)
            elseif(ymd.gt.ny) then
              uxTop = (2*(h1(x,y,z, 8)-h1(x,y,z,13))/(Zliq-1./3.))
              uyTop = (2*(h1(x,y,z,12)-h1(x,y,z,18))/(Zliq-1./3.))
              h2(x,y,z,i) = h1(x,y,z,ibar) - 6*wi(ibar)* (dirx(ibar)*uxTop + diry(ibar)*uyTop)
            elseif(ymd.lt.1) then
              h2(x,y,z,i) = h1(x,y,z,ibar) - 6*wi(ibar)*dirx(ibar)*uxWall
              FxLiqOnWall = FxLiqOnWall + (h2(x,y,z,i) + h1(x,y,z,ibar))*dirx(i)
            elseif((xmd-x_sphere)**2 + (ymd-y_sphere)**2 + (zmd-z_sphere)**2.lt.rad_sphere**2) then
              h2(x,y,z,i) = h1(x,y,z,ibar) - 6*wi(ibar)*dirx(ibar)*uxWall
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
  real*8 :: mass_sphere_SI
  TanStress = (FxLiqOnWall+FxLiqOnSphe)/(nx*nz)
  ZL = TanStress / uxWall
  Zq = (rho0/rho_Liq_SI)*(nu/nu_SI)*(delta_f0_SI/delta_f0)*Zq_SI
  dfc = f0_SI*complex(0,1)/PI *ZL/Zq
  dfc = -conjg(dfc)
  mass_sphere_SI = rho_sphere_SI * 4.*Pi/3. * rad_sphere_SI**3.
  dfc = dfc - 2*ovtorder*f0_SI**2 / Zq_SI * mass_sphere_SI

end subroutine calc_dfc
