module declare
  implicit none
  integer :: i_calc_type
  logical :: B_Save
  integer, parameter :: j_elastic = 1, j_BouncB = 2
  integer, parameter :: ndir = 18, nx_mx = 100, ny_mx = 200, nz_mx = 100, novt_mx = 1 
  integer, parameter :: nPo_mx = nx_mx * ny_mx * nz_mx !, nLi_mx = nx_mx * ny_mx * nz_mx * (ndir + 1)
  integer, parameter :: nPoTop_mx = nx_mx * nz_mx, nPoWal_mx = nPoTop_mx
  integer, parameter :: nPoBor_mx = nx_mx * nz_mx * 3,nPoSph_mx = int(nPo_mx / 5) 
  
  integer, parameter :: dirx(0:18) = (/0,+1,-1, 0, 0, 0, 0,+1,-1,+1,-1, 0, 0,+1,-1,+1,-1, 0, 0/)
  integer, parameter :: dirz(0:18) = (/0, 0, 0, 0, 0,+1,-1, 0, 0,+1,-1,+1,-1, 0, 0,-1,+1,-1,+1/)
  integer, parameter :: diry(0:18) = (/0, 0, 0,+1,-1, 0, 0,+1,-1, 0, 0,+1,-1,-1,+1, 0, 0,+1,-1/)
  integer, parameter :: no_slip(0:18) = (/0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17/)
  real*8,parameter :: wi(0:18) = (/ 1./3.,1./18,1./18.,1./18.,1./18.,1./18.,1./18., &
                           1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36./)
  real*8,parameter :: pi = 4.D0*DATAN(1.D0), rho0 = 1., nu = 1./6., tauinv_Liq = 2./(6.*nu + 1.), omtauinv_Liq = 1. - tauinv_Liq
 
  real*8,parameter :: Zq_SI = 8.8e6, nu_SI = 1e-6, rho0_SI = 1e3, f0_SI = 5e6 
  real*8,parameter :: delta_f0_SI = sqrt(2 * nu_SI/(2*pi*f0_SI));
  real*8           :: f_SI, cell_height_SI, gridreso_SI, delta_f0,f0
  real*8           :: rad_Sph_SI, rho_Sph_SI, G_Sph_SI, eta_Sph_SI, invtandel_to_Gpp_liq_f0
  real*8           :: x_Sph, y_Sph, z_Sph, rad_Sph, G_Sph, nufac_Sph, K_Sph
  real*8           :: sizefac, ncycles, delta, omega, omega2by2 
  real*8           :: tauinv_Sph, omtauinv_Sph, tauinv_Liq_m_omega2by2, tauinv_Sph_m_omega2by2
  complex*16 :: omtauinv_Sph_p_Iomega,Iomega
  
  complex*16 ::               Fx(1:nPo_mx), Fy(1:nPo_mx), Fz(1:nPo_mx)
  complex*16, allocatable :: dr(:),ux(:),uy(:),uz(:),dx(:),dy(:),dz(:),h1(:), h2(:)

  complex*16 :: uxTop(1:nPoTop_mx), uzTop(1:nPoTop_mx)
  complex*16 :: Rimux(1:nPoBor_mx),Rimuy(1:nPoBor_mx),Rimuz(1:nPoBor_mx)
  complex*16 :: FxBor(1:nPoBor_mx,0:ndir),FyBor(1:nPoBor_mx,0:ndir),FzBor(1:nPoBor_mx,0:ndir)
  complex*16 :: Fx3PLSph(1:nPoBor_mx,0:ndir),Fy3PLSph(1:nPoBor_mx,0:ndir),Fz3PLSph(1:nPoBor_mx,0:ndir)
  complex*16 :: Fx3PLWal(1:nPoBor_mx,0:ndir),Fy3PLWal(1:nPoBor_mx,0:ndir),Fz3PLWal(1:nPoBor_mx,0:ndir)
  complex*16 :: FxWalLiq(1:nPoWal_mx,1:5),FxWalSph(1:nPoWal_mx,1:5)
  complex*16 :: FxLiqOnWal,FxLiqOnSph,FxSphOnWal,dfc,ZL,Zq,Zliq
  complex*16 :: dfc_ref(1:novt_mx),dfc_fin(1:novt_mx) 
    
  integer :: Poij(1:nx_mx,1:ny_mx,1:nz_mx)
  integer :: Rimj(1:nx_mx,1:ny_mx,1:nz_mx)
  integer :: Linj(1:nx_mx,1:ny_mx,1:nz_mx,0:ndir)
  integer :: Topj(1:nx_mx,1:nz_mx)
  integer :: nsteps,ovtorder,novt
  integer :: nx,ny,nz
  integer :: nPo,nPoLiq,nPoSph,nPoTop,nPoWalLiq,nPoWalSph,nPoBor,nPoRim,nPo3PL
  integer :: nLi,nLiBor(1:nPoBor_mx),nLi3PLSph(1:nPoBor_mx),nLi3PLWal(1:nPoBor_mx)
  integer :: PoLiqx(1:nPo_mx)   ,PoLiqy(1:nPo_mx   ),PoLiqz(1:nPo_mx)
  integer :: PoSphx(1:nPoSph_mx),PoSphy(1:nPoSph_mx),PoSphz(1:nPoSph_mx)
  integer :: PoTopx(1:nPoTop_mx),PoTopy(1:nPoTop_mx),PoTopz(1:nPoTop_mx)
  integer :: PoBorx(1:nPoBor_mx),PoBory(1:nPoBor_mx),PoBorz(1:nPoBor_mx)
  integer :: Po3PLx(1:nPoBor_mx),Po3PLy(1:nPoBor_mx),Po3PLz(1:nPoBor_mx)
  integer :: PoRimx(1:nPoBor_mx),PoRimy(1:nPoBor_mx),PoRimz(1:nPoBor_mx)
  integer :: PoWalLiqx(1:nPoWal_mx),PoWalLiqy(1:nPoWal_mx),PoWalLiqz(1:nPoWal_mx)
  integer :: PoWalSphx(1:nPoWal_mx),PoWalSphy(1:nPoWal_mx),PoWalSphz(1:nPoWal_mx)
  integer :: permu_Bor(1:ndir+1,1:nPoBor_mx),permu_3PL(1:ndir+1,1:nPoBor_mx)
  integer :: beginning,rate,ending
  integer, parameter :: permu_top_Cov(1:5 ) = (/        4,      8,        13,12,            18/)
  integer, parameter :: permu_top_Liq(1:14) = (/0,1,2,3,  5,6,7,  9,10,11,      14,15,16,17   /)
  integer, parameter :: permu_Wal_Wal(1:5 ) = (/      3,      7,       14,      11,      17   /)
  integer, parameter :: permu_Wal_Liq(1:14) = (/0,1,2,  4,5,6,  8,9,10,   12,13,   15,16,   18/)
  integer :: it,error
  
  complex*16 :: test,Sauerbrey_Contr
  real*8 :: mass_Sph
end module declare

program FT_LBM
  use declare; implicit none
  integer :: iovt,x,y,z
  open(unit = 20,file='freq_ring_in.asc1',action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='fin_freq.asc1',    action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='ux_vs_y.asc1',     action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='uxuy_3D.asc1',     action='write',status='unknown'); close(unit = 20,status='delete');   
  B_Save = .true.
  i_calc_type = j_BouncB; 
  novt = 1; 
  ncycles = 12; 
  sizefac = 3.;
  gridreso_SI = 0.8e-9; 
  rad_Sph_SI = 4e-9; 
  cell_height_SI = 8e-9;
  G_Sph_SI = 1e7; eta_Sph_SI = 1e-3; rho_Sph_SI = 1e3
  ny = int(cell_height_SI/gridreso_SI)
  delta_f0 = delta_f0_SI/gridreso_SI

  f0 = 2.*nu/delta_f0**2 / (2. * pi)
  invtandel_to_Gpp_liq_f0 = G_Sph_SI / (2. * pi * f0_SI * nu_SI * rho0_SI)
  G_Sph = invtandel_to_Gpp_liq_f0 * (2 * pi * f0 * nu * rho0)
  nufac_Sph = (eta_Sph_SI / rho_Sph_SI) / nu_SI
  tauinv_Sph = 2./(6.*nu*nufac_Sph + 1.)
  omtauinv_Sph = 1. - tauinv_sph

  rad_Sph = -1; nx = 1; nz = 1; nPoTop = nx * nz; 
  nLi = nx * ny * nz * (ndir+1); allocate(h1(nLi)); allocate(h2(nLi))
  nPo = nx*ny*nz; allocate(dr(nPo)); 
  allocate(ux(nPo)); allocate(uy(nPo)); allocate(uz(nPo))
  call make_Lists
  do iovt = 1,novt
    call Set_om_dependent_variables(iovt); call init_u0; 
    do it = 1,nsteps
      call compute_VelTop; 
      call Stream_and_Collide; 
      call calc_dfc; 
      if(mod(it,int(nsteps/10)) == 1) then 
         print *,'Liq',it,nsteps,real(dfc),aimag(dfc)
      endif  
    enddo
    dfc_ref(iovt) = dfc;
    open(unit = 20,file='ux_vs_y.asc1', action='write',status='unknown', position='append'); 
    do y = 1,ny 
      write (20,*) y*gridreso_SI*1e9,real(ux(Poij(1,y,1))), aimag(ux(Poij(1,y,1))); 
    enddo  
    close(unit = 20)
  enddo
  deallocate(h1); deallocate(h2); 
  deallocate(dr); deallocate(ux); deallocate(uy); deallocate(uz);



  rad_Sph = rad_Sph_SI * delta_f0 / delta_f0_SI; 
  nx = int(2. * rad_Sph * sizefac); 
  if (mod(nx,2).eq.0) nx = nx + 1;
  nz = nx; x_Sph = nx/2.; z_Sph = nz/2.; y_Sph = 0.
  nPoTop = nx * nz;
  nLi = nx * ny * nz * (ndir+1); allocate(h1(nLi)); allocate(h2(nLi))
  nPo = nx*ny*nz; allocate(dr(nPo)); allocate(ux(nPo)); allocate(uy(nPo)); allocate(uz(nPo))
!  rad_Sph = 1.5 !!!!!!!!!!!!!!!!!!!!
  call system_clock(beginning,rate)
  call make_Lists
  do iovt = 1,novt
    call Set_om_dependent_variables(iovt); 
    call init_u0; 
    do it = 1,nsteps
      call compute_VelTop; 
      if(i_calc_type.eq.j_BouncB) call sphere_dynamics
      call Stream_and_Collide; 
      call calc_dfc; 
      if(mod(it,int(nsteps/10)) == 1) then 
        print *,it,nsteps,real(dfc- dfc_ref(iovt)),aimag(dfc- dfc_ref(iovt))
      endif  
      if (mod(it,int(nsteps/100)) == 1) then 
        open(unit = 20,file='freq_ring_in.asc1',action='write',status='unknown', position='append'); 
        write (20,*)   it,nsteps,real(dfc- dfc_ref(iovt)),aimag(dfc- dfc_ref(iovt)); 
        close(unit = 20);   
      endif  
    enddo
    dfc_fin(iovt) = dfc - dfc_ref(iovt)
  enddo
  deallocate(h1); deallocate(h2);
  open(unit = 20,file='fin_freq.asc1',    action='write',status='unknown', position='append'); 
  write (20,*) cell_height_SI,ovtorder,real(dfc_fin)/ovtorder,aimag(dfc_fin)/ovtorder; close(unit = 20)
  
  open(unit = 20,file='uxuy_3D.asc1',action='write')
  write(20,"(I5,I5,I5)")nx,ny,nz
  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        write(20,"(I5,I5,I5,F12.8,F12.8,F12.8,F12.8,F12.8,F12.8)") x,y,z, &
                                   real(ux(Poij(x,y,z))),aimag(ux(Poij(x,y,z))),&
                                   real(uy(Poij(x,y,z))),aimag(uy(Poij(x,y,z))),&
                                   real(uz(Poij(x,y,z))),aimag(uz(Poij(x,y,z)))
      enddo
    enddo
  enddo
  close(unit = 20)
  deallocate(dr); deallocate(ux); deallocate(uy); deallocate(uz);
end program FT_LBM

logical function in_Sph(x,y,z) 
  use declare; implicit none
  integer :: x,y,z
  in_Sph = ( (y.ge.1) .and. ((x-x_Sph)**2 + (y-y_Sph)**2 + (z-z_Sph)**2 .lt. rad_Sph**2) )
end function in_Sph

logical function in_Liq(x,y,z) 
  use declare; implicit none
  integer :: x,y,z,xmd,ymd,zmd,i
  logical in_Sph,B_mds_not_in_Sph
  B_mds_not_in_Sph = .true.
  do i = 0,ndir
    xmd = mod(nx+x-dirx(i)-1,nx)+1
    ymd = y-diry(i)
    zmd = mod(nz+z-dirz(i)-1,nz)+1
    if ( (in_Sph(xmd,ymd,zmd)) .and. (ymd.ge.1)) then 
      B_mds_not_in_Sph = .false.
    endif  
  enddo  
  in_Liq = ((.not. in_Sph(x,y,z)) .and. (y.gt.1) .and. (B_mds_not_in_Sph))
end function in_Liq

logical function on_WalLiq(x,y,z) 
  use declare; implicit none
  integer :: x,y,z,xmd,ymd,zmd,i
  logical in_Sph,B_mds_not_in_Sph
  B_mds_not_in_Sph = .true.
  do i = 0,ndir
    xmd = mod(nx+x-dirx(i)-1,nx)+1
    ymd = y-diry(i)
    zmd = mod(nz+z-dirz(i)-1,nz)+1
    if ( (in_Sph(xmd,ymd,zmd)) .and. (ymd.ge.1)) then 
      B_mds_not_in_Sph = .false.
    endif  
  enddo  
  on_WalLiq = ((.not. in_Sph(x,y,z)) .and. (y.eq.1) .and. (B_mds_not_in_Sph))
end function on_WalLiq

logical function on_Bor(x,y,z) 
  use declare; implicit none
  integer :: x,y,z,xmd,ymd,zmd,i
  logical in_Sph,B_Bor
  B_Bor = .false.
  do i = 0,ndir
    xmd = mod(nx+x-dirx(i)-1,nx)+1
    ymd = y-diry(i)
    zmd = mod(nz+z-dirz(i)-1,nz)+1
    if ( (.not. in_Sph(x,y,z)) .and. in_Sph(xmd,ymd,zmd) .and. (y.gt.1) ) then 
      B_Bor = .true.
    endif  
  enddo  
  on_Bor = B_Bor
end function on_Bor

logical function on_3PL(x,y,z) 
  use declare; implicit none
  integer :: x,y,z,xmd,ymd,zmd,i
  logical in_Sph,B_3PL
  B_3PL = .false.
  do i = 0,ndir
    xmd = mod(nx+x-dirx(i)-1,nx)+1
    ymd = y-diry(i)
    zmd = mod(nz+z-dirz(i)-1,nz)+1
    if ( (.not. in_Sph(x,y,z)) .and. in_Sph(xmd,ymd,zmd) .and. (ymd.ge.1) .and. (y.eq.1) ) then 
      B_3PL = .true.
    endif  
  enddo  
  on_3PL = B_3PL
end function on_3PL

logical function on_Rim(x,y,z) 
  use declare; implicit none
  integer :: x,y,z,xmd,ymd,zmd,i
  logical in_Sph,B_Rim
  B_Rim = .false.
  do i = 0,ndir
    xmd = mod(nx+x-dirx(i)-1,nx)+1
    ymd = y-diry(i)
    zmd = mod(nz+z-dirz(i)-1,nz)+1
    if ((.not. in_Sph(xmd,ymd,ymd)) .and. in_Sph(x,y,z)) B_Rim = .true.
  enddo  
  on_Rim = B_Rim
end function on_Rim

subroutine Set_om_dependent_variables(iovt)
  use declare; implicit none
  integer :: iovt
  ovtorder = 2 * iovt - 1; f_SI = ovtorder* f0_SI; 
  delta = delta_f0 / sqrt(float(ovtorder)); 
  omega = 2.*nu/delta**2.; omega2by2 = omega**2/2; Zliq = -conjg(sqrt(rho0*complex(0,1)*omega*rho0*nu))
  tauinv_Liq_m_omega2by2 = tauinv_Liq - omega2by2
  tauinv_Sph_m_omega2by2 = tauinv_Sph - omega2by2
  omtauinv_Sph_p_Iomega = omtauinv_Sph + complex(0.,1.)*omega
  Iomega = complex(0.,1.)*omega
  K_Sph = nu * rho0 * omega * 1e4;
  nsteps = int(ncycles / omega * (ny / delta)**2)
end subroutine Set_om_dependent_variables

subroutine make_Lists
  use declare; implicit none
  integer :: x,y,z,xmd,ymd,zmd,i,i_next
  logical :: in_Sph,in_Liq,on_WalLiq,on_Bor,on_3PL,on_Rim
  nPo = 0; nLi = 0; nPoLiq = 0; nPoSph = 0; nPoTop = 0; nPoWalLiq = 0; nPoWalSph = 0; 
  nPoBor = 0; nPo3PL = 0; nPoRim = 0;  
  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        nPo = nPo + 1; 
        Poij(x,y,z) = nPo
        do i = 0,ndir
          nLi = nLi + 1 
          Linj(x,y,z,i) = nLi
        enddo

        if ( (in_Liq(x,y,z)) .and. (y.gt.1) .and. (y.lt.ny) ) then
          nPoLiq = nPoLiq + 1        
          PoLiqx(nPoLiq) = x
          PoLiqy(nPoLiq) = y
          PoLiqz(nPoLiq) = z
        endif   

        if ( on_WalLiq(x,y,z) ) then 
          nPoWalLiq = nPoWalLiq + 1
          PoWalLiqx(nPoWalLiq) = x
          PoWalLiqy(nPoWalLiq) = y
          PoWalLiqz(nPoWalLiq) = z
        endif  

        if ( (y.eq.1) .and. in_Sph(x,y,z) ) then 
          nPoWalSph = nPoWalSph + 1
          PoWalSphx(nPoWalSph) = x
          PoWalSphy(nPoWalSph) = y
          PoWalSphz(nPoWalSph) = z
        endif  

        if (in_Sph(x,y,z)) then
          nPoSph = nPoSph + 1
          PoSphx(nPoSph) = x
          PoSphy(nPoSph) = y
          PoSphz(nPoSph) = z
        endif  

        if (y.eq.ny) then 
          nPoTop = nPoTop + 1
          Topj(x,z) = nPoTop
          PoTopx(nPoTop) = x
          PoTopy(nPoTop) = y
          PoTopz(nPoTop) = z
        endif  

        if (on_Bor(x,y,z)) then
          nPoBor = nPoBor + 1
          PoBorx(nPoBor) = x
          PoBory(nPoBor) = y
          PoBorz(nPoBor) = z
          nLiBor(nPoBor) = 0  
          do i = 0,ndir
            xmd = mod(nx+x-dirx(i)-1,nx)+1
            ymd = y-diry(i)
            zmd = mod(nz+z-dirz(i)-1,nz)+1
            if (in_Sph(xmd,ymd,zmd)) then 
              nLiBor(nPoBor) = nLiBor(nPoBor) + 1
              permu_Bor(nLiBor(nPoBor),nPoBor) = i
            endif
          enddo
          i_next = nLiBor(nPoBor)
          do i = 0,ndir
            xmd = mod(nx+x-dirx(i)-1,nx)+1
            ymd = y-diry(i)
            zmd = mod(nz+z-dirz(i)-1,nz)+1
            if ((ymd.ge.1) .and. (.not.in_Sph(xmd,ymd,zmd)) ) then 
              i_next = i_next + 1
              permu_Bor(i_next,nPoBor) = i
            endif
          enddo
        endif  

        if (on_3PL(x,y,z)) then
          nPo3PL = nPo3PL + 1
          Po3PLx(nPo3PL) = x
          Po3PLy(nPo3PL) = y
          Po3PLz(nPo3PL) = z
          nLi3PLSPh(nPo3PL) = 0  
          nLi3PLWal(nPo3PL) = 0  
          do i = 0,ndir
            xmd = mod(nx+x-dirx(i)-1,nx)+1
            ymd = y-diry(i)
            zmd = mod(nz+z-dirz(i)-1,nz)+1
            if ( (in_Sph(xmd,ymd,zmd)) .and. (ymd.ge.1) ) then 
              nLi3PLSPh(nPo3PL) = nLi3PLSPh(nPo3PL) + 1
              permu_3PL(nLi3PLSPh(nPo3PL),nPo3PL) = i
            endif
          enddo
          i_next = nLi3PLSPh(nPo3PL)
          do i = 0,ndir
            xmd = mod(nx+x-dirx(i)-1,nx)+1
            ymd = y-diry(i)
            zmd = mod(nz+z-dirz(i)-1,nz)+1
            if (ymd.lt.1) then 
              i_next = i_next + 1
              nLi3PLWal(nPo3PL) = nLi3PLWal(nPo3PL) + 1
              permu_3PL(i_next,nPo3PL) = i
            endif
          enddo
          do i = 0,ndir
            xmd = mod(nx+x-dirx(i)-1,nx)+1
            ymd = y-diry(i)
            zmd = mod(nz+z-dirz(i)-1,nz)+1
            if ( ( .not. in_Sph(xmd,ymd,zmd)) .and. (ymd.ge.1) ) then 
              i_next = i_next + 1
              permu_3PL(i_next,nPo3PL) = i
            endif
          enddo
        endif  

        if (on_Rim(x,y,z)) then
          nPoRim = nPoRim + 1
          Rimj(x,y,z) = nPoRim
          PoRimx(nPoRim) = x
          PoRimy(nPoRim) = y
          PoRimz(nPoRim) = z
        endif
      enddo
    enddo
  enddo
end subroutine make_Lists

subroutine init_u0
  use declare; implicit none
  integer :: x,y,z,iPo,iLi
  do iPo = 1,nPo
    dr(iPo) = 0
    ux(iPo) = 0
    uy(iPo) = 0
    uz(iPo) = 0
  enddo  
  if (i_calc_type .eq. j_BouncB) then
    do iPo = 1,nPoSph
      x = PoSphx(iPo)
      y = PoSphy(iPo)
      z = PoSphz(iPo)
      ux(Poij(x,y,z)) = 1
    enddo  
  endif  
  do iLi = 1,nLi
    h1(iLi) = 0
  enddo
end subroutine init_u0

subroutine compute_VelTop
  use declare; implicit none
  integer :: x,z
  do x = 1,nx
    do z = 1,nz
      uxTop(Topj(x,z)) = -2*( h1(Linj(x,ny,z,7 )) - h1(Linj(x,ny,z,14)) )/(Zliq-1./3.)
      uzTop(Topj(x,z)) = -2*( h1(Linj(x,ny,z,11)) - h1(Linj(x,ny,z,17)) )/(Zliq-1./3.)  ! signs?
    enddo  
  enddo
end subroutine compute_VelTop

subroutine relax_Liq(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
  use declare; implicit none
  integer :: x,y,z,i
  complex*16 :: cidotu,heq,drloc,uxloc,uyloc,uzloc,ht(0:ndir)
  cidotu = dirx(i)*uxloc + diry(i)*uyloc + dirz(i)*uzloc
  heq = wi(i) * (drloc + 3.*cidotu)
  h2(Linj(x,y,z,i)) = tauinv_Liq_m_omega2by2 * heq + Iomega * ht(i)   !omtauinv = 0
end subroutine relax_Liq

subroutine relax_Sph(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
  use declare; implicit none
  integer :: x,y,z,i
  complex*16 :: cidotu,heq,drloc,uxloc,uyloc,uzloc,ht(0:ndir)
  cidotu = dirx(i)*uxloc + diry(i)*uyloc + dirz(i)*uzloc
  heq = wi(i) * (drloc + 3.*cidotu)
  h2(Linj(x,y,z,i)) = omtauinv_Sph_p_Iomega * ht(i) + tauinv_Sph_m_omega2by2 * heq 
end subroutine relax_Sph

subroutine calc_dr_uxyz(x,y,z,ht,drloc,uxloc,uyloc,uzloc)
  use declare; implicit none
  integer :: x,y,z,i
  complex*16 :: drloc,ht(0:ndir),uxloc,uyloc,uzloc
  drloc = 0; uxloc = 0; uyloc = 0; uzloc = 0; 
  do i = 0,ndir
    drloc = drloc + ht(i)
    uxloc = uxloc + ht(i) * dirx(i)
    uyloc = uyloc + ht(i) * diry(i)
    uzloc = uzloc + ht(i) * dirz(i)
  enddo  
  if (B_Save) then
    dr(Poij(x,y,z)) = drloc
    ux(Poij(x,y,z)) = uxloc
    uy(Poij(x,y,z)) = uyloc
    uz(Poij(x,y,z)) = uzloc
  endif
end subroutine calc_dr_uxyz

subroutine Apply_Elastic_Forces(x,y,z)
  use declare; implicit none
  integer :: x,y,z,xm1,ym1,zm1,xp1,yp1,zp1,i
  complex*16 :: dxWall,Si
  complex*16 :: d_xx_x, d_yy_x, d_zz_x, d_xy_y, d_xz_z
  complex*16 :: d_xx_y, d_yy_y, d_zz_y, d_xy_x, d_yz_z
  complex*16 :: d_xx_z, d_yy_z, d_zz_z, d_xz_x, d_yz_y
  dxWall = 1/(complex(0,1)*omega)
  dx(Poij(x,y,z)) = ux(Poij(x,y,z)) / (complex(0,1)*omega)     
  dy(Poij(x,y,z)) = uy(Poij(x,y,z)) / (complex(0,1)*omega)     
  dz(Poij(x,y,z)) = uz(Poij(x,y,z)) / (complex(0,1)*omega)
  xm1 = mod(nx+x-1-1,nx)+1
  xp1 = mod(nx+x+1-1,nx)+1
  zm1 = mod(nz+z-1-1,nz)+1
  zp1 = mod(nz+z+1-1,nz)+1
  ym1 = y-1
  yp1 = y+1

  d_xx_x = dx(Poij(xp1,y  ,z  )) - 2*dx(Poij(x,y,z)) + dx(Poij(xm1,y  ,z  ))
  d_zz_x = dx(Poij(x  ,y  ,zp1)) - 2*dx(Poij(x,y,z)) + dx(Poij(x  ,y  ,zm1))
  d_xx_y = dy(Poij(xp1,y  ,z  )) - 2*dy(Poij(x,y,z)) + dy(Poij(xm1,y  ,z  ))
  d_zz_y = dy(Poij(x  ,y  ,zp1)) - 2*dy(Poij(x,y,z)) + dy(Poij(x  ,y  ,zm1))
  d_xx_z = dz(Poij(xp1,y  ,z  )) - 2*dz(Poij(x,y,z)) + dz(Poij(xm1,y  ,z  ))
  d_zz_z = dz(Poij(x  ,y  ,zp1)) - 2*dz(Poij(x,y,z)) + dz(Poij(x  ,y  ,zm1))
  d_xz_z = dz(Poij(xp1,y  ,zp1)) - dz(Poij(xp1,y  ,zm1)) - dz(Poij(xm1,y  ,zp1)) + dz(Poij(xm1,y  ,zm1))
  d_xz_x = dx(Poij(xp1,y  ,zp1)) - dx(Poij(xp1,y  ,zm1)) - dx(Poij(xm1,y  ,zp1)) + dx(Poij(xm1,y  ,zm1))

  if (y.gt.1) then 
    d_yy_x = dx(Poij(x  ,yp1,z  )) - 2*dx(Poij(x,y,z)) + dx(Poij(x  ,ym1,z  ))
    d_yy_y = dy(Poij(x  ,yp1,z  )) - 2*dy(Poij(x,y,z)) + dy(Poij(x  ,ym1,z  ))
    d_yy_z = dz(Poij(x  ,yp1,z  )) - 2*dz(Poij(x,y,z)) + dz(Poij(x  ,ym1,z  ))
    d_xy_x = dx(Poij(xp1,yp1,z  )) - dx(Poij(xp1,ym1,z  )) - dx(Poij(xm1,yp1,z  )) + dx(Poij(xm1,ym1,z  ))
    d_xy_y = dy(Poij(xp1,yp1,z  )) - dy(Poij(xp1,ym1,z  )) - dy(Poij(xm1,yp1,z  )) + dy(Poij(xm1,ym1,z  ))
    d_yz_z = dy(Poij(x  ,yp1,zp1)) - dz(Poij(x  ,ym1,zp1)) - dz(Poij(x  ,yp1,zm1)) + dz(Poij(x  ,ym1,zm1))
  endif

  if (y.eq.1) then 
    d_yy_x = dx(Poij(x  ,yp1,z  )) - 2*dx(Poij(x,y,z)) + 0
    d_yy_y = dy(Poij(x  ,yp1,z  )) - 2*dy(Poij(x,y,z)) + 0
    d_yy_z = dz(Poij(x  ,yp1,z  )) - 2*dz(Poij(x,y,z)) + 0
    d_xy_x = dx(Poij(xp1,yp1,z  )) - 0 - dx(Poij(xm1,yp1,z  )) + 0
    d_xy_y = dy(Poij(xp1,yp1,z  )) - 0 - dy(Poij(xm1,yp1,z  )) + 0
    d_yz_z = dy(Poij(x  ,yp1,zp1)) - 0 - dz(Poij(x  ,yp1,zm1)) + 0
  endif  
        
  Fx(Poij(x,y,z)) = (3 * K_Sph - 2 * G_Sph) * (d_xx_x + d_xy_y + d_xz_z) + & 
                    G_Sph * (d_xx_x + d_xy_y + d_xz_z + d_xx_x + d_yy_x + d_zz_x)
  Fy(Poij(x,y,z)) = (3 * K_Sph - 2 * G_Sph) * (d_xy_x + d_yy_y + d_yz_z) + & 
                    G_Sph * (d_xy_x + d_yy_y + d_yz_z + d_yy_y + d_yy_y + d_zz_y)
  Fz(Poij(x,y,z)) = (3 * K_Sph - 2 * G_Sph) * (d_xz_x + d_yz_y + d_zz_z) + &
                    G_Sph * (d_xz_x + d_yz_y + d_zz_z + d_xx_z + d_yy_z + d_zz_z)
  do i = 0,ndir
    Si = (1.-tauinv_Sph/2.)*3.*(dirx(i)*Fx(Poij(x,y,z)) + diry(i)*Fy(Poij(x,y,z)) + dirz(i)*Fz(Poij(x,y,z)))
    h1(Linj(x,y,z,i)) = h1(Linj(x,y,z,i)) + Si
  enddo
end subroutine Apply_Elastic_Forces

subroutine Stream_and_Collide
  use declare; implicit none
  integer :: x,y,z,i,iloc,ibar,xmd,ymd,zmd,iPo,iLi,iRim
  complex*16 :: ht(0:ndir),drloc,uxloc,uyloc,uzloc,cidotu
  do iPo = 1,nPoBor
    do iloc = 0,ndir
      FxBor(iPo,iloc) = 0
      FyBor(iPo,iloc) = 0
      FzBor(iPo,iloc) = 0
    enddo
  enddo  

  do iPo = 1,nPoWalLiq
    x = PoWalLiqx(iPo)
    y = PoWalLiqy(iPo)
    z = PoWalLiqz(iPo)
    do iloc = 1,5
      i = permu_Wal_Wal(iloc)
      ibar = no_slip(i)
      ht(i) = h1(Linj(x,y,z,ibar)) - 6 * wi(ibar) * dirx(ibar)
      FxWalLiq(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * dirx(i)
    enddo  
    do iloc = 1,14
      i = permu_Wal_Liq(iloc)
      xmd = mod(nx+x-dirx(i)-1,nx)+1
      ymd = y-diry(i)
      zmd = mod(nz+z-dirz(i)-1,nz)+1
      ht(i) = h1(Linj(xmd,ymd,zmd,i))
    enddo  
    call calc_dr_uxyz(x,y,z,ht,drloc,uxloc,uyloc,uzloc)
    do i = 0,ndir
      call relax_Liq(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
    enddo  
  enddo  
  
  do iPo = 1,nPoLiq
    x = PoLiqx(iPo)
    y = PoLiqy(iPo)
    z = PoLiqz(iPo)
    do i = 0,ndir
      xmd = mod(nx+x-dirx(i)-1,nx)+1
      ymd = y-diry(i)
      zmd = mod(nz+z-dirz(i)-1,nz)+1
      ht(i) = h1(Linj(xmd,ymd,zmd,i))
    enddo 
    call calc_dr_uxyz(x,y,z,ht,drloc,uxloc,uyloc,uzloc)
    do i = 0,ndir
      call relax_Liq(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
    enddo  
  enddo  

  if (i_calc_type.eq.j_elastic) then
    do iPo = 1,nPoWalSph
      x = PoWalSphx(iPo)
      y = PoWalSphy(iPo)
      z = PoWalSphz(iPo)
      do iloc = 1,5
        i = permu_Wal_wal(iloc)
        ibar = no_slip(i)
        ht(i) = h1(Linj(x,y,z,ibar)) - 6 * wi(ibar) * dirx(ibar)
        FxWalSph(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * dirx(i)
      enddo  
      do iloc = 1,14
        i = permu_Wal_Liq(iloc)
        xmd = mod(nx+x-dirx(i)-1,nx)+1
        ymd = y-diry(i)
        zmd = mod(nz+z-dirz(i)-1,nz)+1
        ht(i) = h1(Linj(xmd,ymd,zmd,i))
      enddo  
      call calc_dr_uxyz(x,y,z,ht,drloc,uxloc,uyloc,uzloc)
      do i = 0,ndir
        call relax_Sph(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
      enddo  
      call Apply_Elastic_Forces(x,y,z)
    enddo
  endif  

  if (i_calc_type.eq.j_elastic) then
    do iPo = 1,nPoSph
      x = PoSphx(iPo)
      y = PoSphy(iPo)
      z = PoSphz(iPo)
      do i = 0,ndir
        xmd = mod(nx+x-dirx(i)-1,nx)+1
        ymd = y-diry(i)
        zmd = mod(nz+z-dirz(i)-1,nz)+1
        ht(i) = h1(Linj(xmd,ymd,zmd,i))
      enddo  
      call calc_dr_uxyz(x,y,z,ht,drloc,uxloc,uyloc,uzloc)
      do i = 0,ndir
        call relax_Sph(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
      enddo  
      call Apply_Elastic_Forces(x,y,z)
    enddo  
  endif  

  do iPo = 1,nPoBor
    x = PoBorx(iPo)
    y = PoBory(iPo)
    z = PoBorz(iPo)
    do iloc = 1,nLiBor(iPo)
      i = permu_Bor(iloc,iPo)
      xmd = mod(nx+x-dirx(i)-1,nx)+1
      ymd = y-diry(i)
      zmd = mod(nz+z-dirz(i)-1,nz)+1
      if (i_calc_type.eq.j_BouncB) then
        ibar = no_slip(i)
        iRim = Rimj(xmd,ymd,zmd)
        cidotu = dirx(ibar)*Rimux(iRim) + diry(ibar)*Rimuy(iRim) + diry(ibar)*Rimuz(iRim)
        ht(i) = h1(Linj(x,y,z,ibar)) - 6 * wi(ibar) * cidotu
        FxBor(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * dirx(i)
        FyBor(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * diry(i)
        FzBor(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * dirz(i)
      else   
        ht(i) = h1(Linj(xmd,ymd,zmd,i))
      endif
    enddo  
    do iloc = nLiBor(iPo)+1,ndir+1
      i = permu_Bor(iloc,iPo)
      xmd = mod(nx+x-dirx(i)-1,nx)+1
      ymd = y-diry(i)
      zmd = mod(nz+z-dirz(i)-1,nz)+1
      ht(i) = h1(Linj(xmd,ymd,zmd,i))
    enddo  
    call calc_dr_uxyz(x,y,z,ht,drloc,uxloc,uyloc,uzloc)
    do i = 0,ndir
      call relax_Liq(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
    enddo  
  enddo  

  do iPo = 1,nPo3PL
    x = Po3PLx(iPo)
    y = Po3PLy(iPo)
    z = Po3PLz(iPo)
    do iloc = 1,nLi3PLSPh(iPo)
      i = permu_3PL(iloc,iPo)
      xmd = mod(nx+x-dirx(i)-1,nx)+1
      ymd = y-diry(i)
      zmd = mod(nz+z-dirz(i)-1,nz)+1
      if (i_calc_type.eq.j_BouncB) then
        ibar = no_slip(i)
        iRim = Rimj(xmd,ymd,zmd)
        cidotu = dirx(ibar)*Rimux(iRim) + diry(ibar)*Rimuy(iRim) + diry(ibar)*Rimuz(iRim)
        ht(i) = h1(Linj(x,y,z,ibar)) - 6 * wi(ibar) * cidotu
        Fx3PLSph(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * dirx(i)
        Fy3PLSph(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * diry(i)
        Fz3PLSph(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * dirz(i)
      else   
        ht(i) = h1(Linj(xmd,ymd,zmd,i))
      endif
    enddo  
    do iloc = nLi3PLSPh(iPo) + 1, nLi3PLSPh(iPo) + nLi3PLWal(iPo)
      i = permu_3PL(iloc,iPo)
      xmd = mod(nx+x-dirx(i)-1,nx)+1
      ymd = y-diry(i)
      zmd = mod(nz+z-dirz(i)-1,nz)+1
      if (i_calc_type.eq.j_BouncB) then
        ibar = no_slip(i)
        ht(i) = h1(Linj(x,y,z,ibar)) - 6 * wi(ibar) * dirx(ibar)
        Fx3PLWal(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * dirx(i)
        Fy3PLWal(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * diry(i)
        Fz3PLWal(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * dirz(i)
      else   
        ht(i) = h1(Linj(xmd,ymd,zmd,i))
      endif
    enddo  
    do iloc = nLi3PLSPh(iPo) + nLi3PLWal(iPo) + 1,ndir+1
      i = permu_3PL(iloc,iPo)
      xmd = mod(nx+x-dirx(i)-1,nx)+1
      ymd = y-diry(i)
      zmd = mod(nz+z-dirz(i)-1,nz)+1
      ht(i) = h1(Linj(xmd,ymd,zmd,i))
    enddo  
    call calc_dr_uxyz(x,y,z,ht,drloc,uxloc,uyloc,uzloc)
    do i = 0,ndir
      call relax_Liq(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
    enddo  
  enddo  

  do iPo = 1,nPoTop
    x = PoTopx(iPo)
    y = PoTopy(iPo)
    z = PoLiqz(iPo)
    do iloc = 1,5
      i = permu_top_cov(iloc)
      ibar = no_slip(i)
      cidotu = dirx(ibar)*uxTop(Topj(x,z)) + dirz(ibar)*uzTop(Topj(x,z))
      ht(i) = h1(Linj(x,y,z,ibar)) - 6 * wi(ibar) * cidotu
    enddo  
    do iloc = 1,14
      i = permu_top_liq(iloc)
      xmd = mod(nx+x-dirx(i)-1,nx)+1
      ymd = y-diry(i)
      zmd = mod(nz+z-dirz(i)-1,nz)+1
      ht(i) = h1(Linj(xmd,ymd,zmd,i))
    enddo  
    call calc_dr_uxyz(x,y,z,ht,drloc,uxloc,uyloc,uzloc)
    do i = 0,ndir
      call relax_Liq(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
    enddo  
  enddo

  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        do i = 0,ndir
          h1(Linj(x,y,z,i)) = h2(Linj(x,y,z,i))
        enddo
      enddo
    enddo
  enddo  

!  do iLi = 1,nLi
!    h1(iLi) = h2(iLi)
!  enddo
end subroutine Stream_and_Collide

subroutine sphere_dynamics
  use declare; implicit none
  integer :: iPo,iloc
  mass_Sph = (rho_Sph_SI/rho0_SI) * rho0 * 1./2. * 4.*Pi/3. * rad_Sph**3. !truncation
  FxLiqOnSph = 0.
  do iPo = 1,nPoBor
    do iloc = 1,nLiBor(iPo)
      FxLiqOnSph = FxLiqOnSph + FxBor(iPo,iloc)
    enddo
  enddo
  do iPo = 1,nPo3PL
    do iloc = 1,nLi3PLSph(iPo)
      FxLiqOnSph = FxLiqOnSph + FxBor(iPo,iloc)
    enddo
  enddo
  FxSphOnWal = FxLiqOnSph - complex(0.,1.) * omega * mass_Sph 
  do iPo = 1,nPoRim
    Rimux(iPo) = 1;
    Rimuy(iPo) = 0;
    Rimuz(iPo) = 0;
  enddo  
end subroutine sphere_dynamics

subroutine calc_dfc
  use declare; implicit none
  integer :: iPo,iloc
  FxLiqOnWal = 0.
  do iPo = 1,nPoWalLiq
    do iloc = 1,5
      FxLiqOnWal = FxLiqOnWal + FxWalLiq(iPo,iloc)     
    enddo
  enddo
    do iPo = 1,nPo3PL
    do iloc = nLi3PLSph(iPo) + 1, nLi3PLSph(iPo) + nLi3PLWal(iPo)
      FxLiqOnWal = FxLiqOnWal + Fx3PLWal(iPo,iloc)
    enddo
  enddo

  if (i_calc_type.eq.j_elastic) then
    do iPo = 1,nPoWalSph
      do iloc = 1,5
        FxLiqOnWal = FxLiqOnWal + FxWalSph(iPo,iloc)     
      enddo
    enddo  
  else
    FxLiqOnWal = FxLiqOnWal + FxSphOnWal
  endif  
  ZL = (FxLiqOnWal)/(nx*nz)
  Zq = (rho0/rho0_SI)*(nu/nu_SI)*(delta_f0_SI/delta_f0)*Zq_SI
  Sauerbrey_Contr = complex(0,1) * FxSphOnWal/(nx*nz) * f0_SI /(pi*Zq)
  Sauerbrey_Contr = -conjg(Sauerbrey_Contr)
  dfc = f0_SI * complex(0,1)/(pi*Zq) * ZL
  dfc = -conjg(dfc)
end subroutine calc_dfc


