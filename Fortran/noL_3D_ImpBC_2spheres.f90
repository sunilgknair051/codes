module declare
  implicit none
  integer, parameter :: ndir = 18, nx_mx1 = 70, nz_mx1 = 70
  integer, parameter :: ny_mx1 = 50, nx_mx2 = 70, nz_mx2 = 70, ny_mx2 = 50, novt_mx = 5
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
  real*8           :: rad_sphere_SI1, cell_height_SI1, gridreso_SI1,rad_sphere_SI2, cell_height_SI2, gridreso_SI2,delta_f0

  integer :: beginning, rate, ending,nx1,ny1,nz1,nx2,ny2,nz2,nsize,iovt,novt,nsteps,ovtorder,ovtordmin,it
  real*8 :: x_sphere1,y_sphere1,z_sphere1,Sizefac,delta,omega,omega2by2,f_SI,ncycles, rad_sphere1, rad_sphere2
  complex*16 :: dr1(1:nx_mx1,1:ny_mx1,1:nz_mx1)
  complex*16 :: ux1(1:nx_mx1,1:ny_mx1,1:nz_mx1)
  complex*16 :: uy1(1:nx_mx1,1:ny_mx1,1:nz_mx1)
  complex*16 :: uz1(1:nx_mx1,1:ny_mx1,1:nz_mx1)
  !complex*16 :: dr2(1:nx_mx2,1:ny_mx2,1:nz_mx2), ux2(1:nx_mx2,1:ny_mx2,1:nz_mx2), uy2(1:nx_mx2,1:ny_mx2,1:nz_mx2), uz2(1:nx_mx2,1:ny_mx2,1:nz_mx2)
  complex*16 :: h1(1:nx_mx1,1:ny_mx1,1:nz_mx1,0:ndir), h2(1:nx_mx1,1:ny_mx1,1:nz_mx1,0:ndir)
  complex*16 :: FxLiqOnSphe, FxLiqOnWall, dfc,dfc_upb,TanStress,ZL,Zq,uxTop,uyTop,Zliq
  complex*16 :: dfc_ref(1:novt_mx),dfc_fin(1:novt_mx),dfc_fin_upb(1:novt_mx)
  complex*16 :: uxtop_ref(1:novt_mx)
end module declare

program FT_LBM
  use declare; implicit none
  integer :: x,y,z,isize
  open(unit = 20,file='complex_freq.asc1', action='write',status='unknown'); close(unit = 20,status='delete');
  open(unit = 21,file='fin_freq.asc1',     action='write',status='unknown'); close(unit = 21,status='delete');
  open(unit = 22,file='fin_freq_upb.asc1', action='write',status='unknown'); close(unit = 22,status='delete');


  nsize = 10; novt = 5; !nsteps = 10001 !increased overtone orders
  ncycles = 0.2;  ovtordmin = 7
  rad_sphere_SI1 = 5e-9; cell_height_SI1= 8e-9; gridreso_SI1 = 1e-9; !REMEMBER!! values increase once the ratio of rad/cellheight is not 1:4
  rad_sphere_SI2 = 5e-9; cell_height_SI2 = 8e-9; gridreso_SI2 = 1e-9;
  ny = int(cell_height_SI1/gridreso_SI1)
  delta_f0 = delta_f0_SI * float(ny)/cell_height_SI1
  do iovt = 1,novt
    nx = 1; nz = 1; rad_sphere1 = 0;
    ovtorder = 2 * iovt - 1 + (ovtordmin - 1); f_SI = ovtorder* f0_SI;
    delta = delta_f0 / sqrt(float(ovtorder));
    omega = 2.*nu/delta**2.; omega2by2 = omega**2/2; Zliq = -conjg(sqrt(rho0*complex(0,1)*omega*rho0*nu))
    nsteps = int(ncycles*2*PI/omega)
    call init_u0;
    do it = 1,nsteps
      call LB_Stream;
      call compute_dr_u; call collide; call calc_dfc;
!      if ((it.ge.1500).and.(mod(it,100)==0)) then
      if ((it.ge.500).and.(mod(it,10) == 1)) then
        open(unit = 20,file='complex_freq.asc1',action='write',status='unknown', position='append');
        write (20,*) isize,nsize,sizefac,ovtorder,it,real(dfc),aimag(dfc),real(dfc_upb),aimag(dfc_upb);close(unit = 20);
      endif
      if(mod(it,2000) == 1) print *,it,nsteps,ovtorder,real(dfc),aimag(dfc)
    enddo
    dfc_ref(iovt) = dfc
    uxtop_ref(iovt) = uxtop
  enddo

  rad_sphere1 = rad_sphere_SI1 * delta_f0 / delta_f0_SI;
  rad_sphere2 = rad_sphere_SI2 * delta_f0 / delta_f0_SI;
  print *,'rad_sphere , delta, ny',rad_sphere1,delta,ny
  do isize = 1,nsize
    call system_clock(beginning,rate)
    if (nsize.gt.1) Sizefac = 1.5 + (isize-1.) / (nsize-1.)*2.
    if (nsize.eq.1) Sizefac = 1.5
    nx1 = int(2 * rad_sphere1* Sizefac); nz1 = nx1; x_sphere1 = nx1/2.; z_sphere1 = nz1/2.
    nx2= int(2 * rad_sphere2* Sizefac); nz2 = nx2; x_sphere2 = nx2/2.; z_sphere2 = nz2/2.
    do iovt = 1,novt
      ovtorder = 2 * iovt - 1 + (ovtordmin - 1);
      delta = delta_f0 / sqrt(float(ovtorder));
      omega = 2.*nu/delta**2.; omega2by2 = omega**2/2; Zliq = -conjg(sqrt(rho0*complex(0,1)*omega*rho0*nu))
      nsteps = int(ncycles*2*PI/omega)
      call init_u0;
      do it = 1,nsteps
        call LB_Stream;
        call compute_dr_u; call collide; call calc_dfc;
        if ((it.ge.500).and.(mod(it,10) == 1)) then
          open(unit = 20,file='complex_freq.asc1',action='write',status='unknown', position='append');
          write (20,*) isize,nsize,sizefac,ovtorder,it,real(dfc),aimag(dfc),real(dfc_upb),aimag(dfc_upb);close(unit = 20);
!          write (20,*) it,real(dfc),aimag(dfc);close(unit = 20);
        endif
        if(mod(it,500) == 1) print *,'bb',it,nsteps,ovtorder,real(dfc-dfc_ref(iovt)),aimag(dfc-dfc_ref(iovt))
        if(mod(it,500) == 1) print *,'upb',it,nsteps,ovtorder,real(dfc_upb),aimag(dfc_upb)
      enddo
      dfc_fin(iovt) = dfc - dfc_ref(iovt)
      dfc_fin_upb(iovt) = dfc_upb
      print *,'iovt',iovt,dfc,dfc_ref(iovt),dfc_fin(iovt)
      print *,'fin',isize,nsize,ovtorder,real(dfc_fin/ovtorder),aimag(dfc_fin/ovtorder)
      call system_clock(ending); print *,"Elapsed time is in seconds --->",real(ending - beginning)/real(rate)
    enddo
    open(unit = 21,file='fin_freq.asc1',    action='write',status='unknown', position='append');
    write (21,*) isize,nsize,pi/4./Sizefac**2,real(dfc_fin)/ovtorder,aimag(dfc_fin)/ovtorder; close(unit = 21)
    open(unit = 22,file='fin_freq_upb.asc1',    action='write',status='unknown', position='append');
    write (22,*) isize,nsize,pi/4./Sizefac**2,real(dfc_fin_upb)/ovtorder,aimag(dfc_fin_upb)/ovtorder; close(unit = 22)
  enddo
  ! open(unit = 30,file='uxuy_3D.asc1',action='write')
  ! write(30,"(I5,I5,I5)")nx,ny,nz
  ! do x = 1,nx
  !   do y = 1,ny
  !     do z = 1,nz
  !       write(30,"(I5,I5,I5,F12.8,F12.8,F12.8,F12.8,F12.8,F12.8)")x,y,z,real(ux(x,y,z)),aimag(ux(x,y,z)),&
  !                                      real(uy(x,y,z)),aimag(uy(x,y,z)),real(uz(x,y,z)),aimag(uz(x,y,z))
  !     enddo
  !   enddo
  ! enddo
  ! close(unit = 30)
end program FT_LBM

subroutine init_u0
  use declare; implicit none
  integer :: xpos1,ypos1,zpos1,i
  do xpos1 = 1,nx1!; xpos2 = 1,nx2
    do ypos1 = 1,ny1!; ypos2 = 1,ny2
      do zpos1 = 1,nz1!; zpos2 = 1,nz2
        do i = 0,ndir
          h1(xpos1,ypos1,zpos1,i) = 0
        end do
      end do
    end do
  end do
end subroutine init_u0

subroutine LB_Stream
  use declare; implicit none
  integer :: xpos1,ypos1,zpos1,xmd,ymd,zmd,i,ibar
  FxLiqOnWall = 0; FxLiqOnSphe = 0;
  do xpos1 = 1,nx1
    do ypos1 = 1,ny1
      do zpos1 = 1,nz1
        do i = 0,ndir
          if ((ypos1.le.ny1).and.(ypos1.ge.1)).and.((xpos1-x_sphere1)**2 +(ypos1-y_sphere1)**2 +(zpos1-z_sphere1)**2.gt.rad_sphere1**2)) then
            ibar = no_slip(i)
            xmd = mod(nx1+xpos1-dirx(i)-1,nx1)
            ymd = ypos1-diry(i)
            zmd = mod(nz1+zpos1-dirz(i)-1,nz1)+1
            if((ymd.le.ny).and.(ymd.ge.1)).and.((xmd-x_sphere1)**2 +(ymd-y_sphere1)**2 +(zmd-z_sphere1)**2.ge.rad_sphere1**2)) then
              h2(xpos1,ypos1,zpos1,i) = h1(xmd,ymd,zmd,i)
            elseif(ymd.gt.ny) then
              uxtop = (2*(h1(xpos1,ypos1,zpos1, 8)-h1(xpos1,ypos1,zpos1,13))/(Zliq-1./3.))
              uytop = (2*(h1(xpos1,ypos1,zpos1,12)-h1(xpos1,ypos1,zpos1,18))/(Zliq-1./3.))
              h2(xpos1,ypos1,zpos1,i) = h1(xpos1,ypos1,zpos1,ibar) - 6*wi(ibar)* (dirx(ibar)*uxtop + diry(ibar)*uytop)
            elseif(ymd.lt.1) then
              h2(xpos1,ypos1,zpos1,i) = h1(xpos1,ypos1,zpos1,ibar) - 6*wi(ibar)*dirx(ibar)*uxWall
              FxLiqOnWall = FxLiqOnWall + (h2(xpos1,ypos1,zpos1,i) + h1(xpos1,ypos1,zpos1,ibar))*dirx(i)
            elseif((xmd-x_sphere1)**2 + (ymd-y_sphere1)**2 + (zmd-z_sphere1)**2.lt.rad_sphere1**2) then
              h2(xpos1,ypos1,zpos1,i) = h1(xpos1,ypos1,zpos1,ibar) - 6*wi(ibar)*dirx(ibar)*uxWall
              FxLiqOnSphe = FxLiqOnSphe + (h2(xpos1,ypos1,zpos1,i) + h1(xpos1,ypos1,zpos1,ibar))*dirx(i)
            endif
          endif
        end do
      end do
    end do
  end do
end subroutine LB_Stream

subroutine compute_dr_u
  use declare; implicit none
  integer :: xpos1,ypos1,zpos1,i
  do xpos1 = 1,nx1
    do ypos1 = 1,ny1
      do zpos1 = 1,nz1
        if ((ypos1.le.ny1).and.(ypos1.ge.1).and.((xpos1-x_sphere1)**2 +(ypos1-y_sphere1)**2 +(zpos1-z_sphere1)**2.ge.rad_sphere1**2)) then
          dr1(xpos1,ypos1,zpos1) = complex(0,0)
          ux1(xpos1,ypos1,zpos1) = complex(0,0)
          uy1(xpos1,ypos1,zpos1) = complex(0,0)
          uz1(xpos1,ypos1,zpos1) = complex(0,0)
          do i = 0,ndir
            dr1(xpos1,ypos1,zpos1) = dr1(xpos1,ypos1,zpos1) + h2(xpos1,ypos1,zpos1,i)
            ux1(xpos1,ypos1,zpos1) = ux1(xpos1,ypos1,zpos1) + h2(xpos1,ypos1,zpos1,i)*dirx(i)
            uy1(xpos1,ypos1,zpos1) = uy1(xpos1,ypos1,zpos1) + h2(xpos1,ypos1,zpos1,i)*diry(i)
            uz1(xpos1,ypos1,zpos1) = uz1(xpos1,ypos1,zpos1) + h2(xpos1,ypos1,zpos1,i)*dirz(i)
          enddo
        endif
      enddo
    enddo
  enddo
end subroutine compute_dr_u

subroutine collide
  use declare; implicit none
  complex*16 :: heq,cidotu
  integer :: xpos1,ypos1,zpos1,i
  do xpos1 = 1,nx1
    do ypos1 = 1,ny1
      do zpos1 = 1,nz1
        if ((ypos1.le.ny).and.(ypos1.ge.1).and.((xpos1-x_sphere1)**2 +(ypos1-y_sphere)**2 +(zpos1-z_sphere)**2.ge.rad_sphere1**2)) then
          do i = 0,ndir
            cidotu = dirx(i)*ux1(xpos1,ypos1,zpos1) + diry(i)*uy1(xpos1,ypos1,zpos1) + dirz(i)*uz1(xpos1,ypos1,zpos1)
            heq = wi(i) * (dr1(xpos1,ypos1,zpos1) + 3*cidotu)
!            h1(x,y,z,i) = omtauinv * h2(x,y,z,i) + tauinv*heq +complex(0,1)*omega*h2(x,y,z,i) - omega**2/2*heq
            h1(xpos1,ypos1,zpos1,i) = tauinv*heq + complex(0,1)*omega*h2(xpos1,ypos1,zpos1,i) - omega2by2 *heq  ! omtauinv is zero
          end do
        endif
      enddo
    enddo
  enddo
end subroutine collide

subroutine calc_dfc
  use declare; implicit none
  real*8 :: mass_sphere_SI1, mass_sphere_SI2
  TanStress = (FxLiqOnWall+FxLiqOnSphe)/(nx1*nz1)
  ZL = TanStress / uxWall
  Zq = (rho0/rho_Liq_SI)*(nu/nu_SI)*(delta_f0_SI/delta_f0)*Zq_SI
  dfc = f0_SI*complex(0,1)/PI *ZL/Zq
  dfc = -conjg(dfc)
  mass_sphere_SI1 = rho_sphere_SI * 4.*Pi/3. * rad_sphere_SI1**3.
  mass_sphere_SI2 = rho_sphere_SI * 4.*Pi/3. * rad_sphere_SI2**3.
  dfc = dfc !- 2*ovtorder*f0_SI**2 / Zq_SI * mass_sphere_SI
  if (abs(uxtop_ref(iovt)).gt.0) then
    dfc_upb = f0_SI*complex(0,1)/(PI * Zq)*Zliq*(uxtop/uxtop_ref(iovt)-1)
  endif

end subroutine calc_dfc
