module declare
  implicit none
  logical :: B_Save
  integer, parameter :: j_BounceB = 1, j_elastic = 2, ndir = 18
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
  integer :: i_calc_type,nsteps_mx,ovtorder,novt
  integer :: nx,ny,nz,nLi,nPo,nPoLiq,nPoSph,nPoTop,nPoWalLiq,nPoWalSph,nPoBor,nPoRim,nPo3PL
  complex*16 :: FxLiqOnWal,FxLiqOnSph,FxSphOnWal,dfc,ZL,Zq,Zliq
  
  complex*16, allocatable :: Fx(:), Fy(:), Fz(:)
  complex*16, allocatable :: dr(:),ux(:),uy(:),uz(:),dx(:),dy(:),dz(:),h1(:), h2(:)

  complex*16, allocatable :: uxTop(:), uzTop(:)
  complex*16, allocatable :: Rimux(:),Rimuy(:),Rimuz(:)
  complex*16, allocatable :: FxBor(:,:),FyBor(:,:),FzBor(:,:)
  complex*16, allocatable :: Fx3PLSph(:,:),Fy3PLSph(:,:),Fz3PLSph(:,:)
  complex*16, allocatable :: Fx3PLWal(:,:),Fy3PLWal(:,:),Fz3PLWal(:,:)
  complex*16, allocatable :: FxWalLiq(:,:),FxWalSph(:,:)
  complex*16, allocatable :: dfc_ref(:),dfc_fin(:) 
    
  integer, allocatable :: Poij(:,:,:),Linj(:,:,:,:),Rimj(:,:,:),Topj(:,:)
  integer, allocatable :: nLiBor(:),nLi3PLSph(:),nLi3PLWal(:)
  integer, allocatable :: PoLiqx(:),PoLiqy(:),PoLiqz(:)
  integer, allocatable :: PoSphx(:),PoSphy(:),PoSphz(:)
  integer, allocatable :: PoTopx(:),PoTopy(:),PoTopz(:)
  integer, allocatable :: PoBorx(:),PoBory(:),PoBorz(:)
  integer, allocatable :: Po3PLx(:),Po3PLy(:),Po3PLz(:)
  integer, allocatable :: PoRimx(:),PoRimy(:),PoRimz(:)
  integer, allocatable :: PoWalLiqx(:),PoWalLiqy(:),PoWalLiqz(:)
  integer, allocatable :: PoWalSphx(:),PoWalSphy(:),PoWalSphz(:)
  integer, allocatable :: permu_Bor(:,:),permu_3PL(:,:)
  integer, allocatable :: nsteps(:)
  integer :: beginning,rate,ending
  integer, parameter :: permu_Top_BBa(1:5 ) = (/        4,      8,        13,12,            18/)
  integer, parameter :: permu_Top_Liq(1:14) = (/0,1,2,3,  5,6,7,  9,10,11,      14,15,16,17   /)
  integer, parameter :: permu_Wal_BBA(1:5 ) = (/      3,      7,       14,      11,      17   /)
  integer, parameter :: permu_Wal_Liq(1:14) = (/0,1,2,  4,5,6,  8,9,10,   12,13,   15,16,   18/)
  integer :: it,error
  
  complex*16 :: test,Sauerbrey_Contr
  real*8 :: mass_Sph,target_acc_Hz
end module declare

program FT_LBM
  use declare; implicit none
  integer :: iovt,x,y,z,iPo
 
  real * 8 :: errorprev
  open(unit = 20,file='freq_ring_in.asc1',action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='fin_freq.asc1',    action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='ux_vs_y.asc1',     action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='uxuy_3D.asc1',     action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 21,file='config.asc1',      action='read',status='unknown'); 
  read(21,*) i_calc_type
  close(unit = 20);   
!  i_calc_type = j_BounceB; 
  print *,'i_calc_type',i_calc_type
  novt = 1; allocate( nsteps(novt) )
  ncycles = 12; 
  sizefac = 3.;
  gridreso_SI = 1.5e-9; 
  rad_Sph_SI = 4e-9; 
  cell_height_SI = 12e-9;
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
  target_acc_Hz = 1; errorprev = 5;

  call Count_and_Allocate
  call make_Lists
  do iovt = 1,novt
    call Set_om_dependent_variables(iovt); 
    nsteps(iovt) = nsteps_mx
    do it = 1,nsteps_mx
      call compute_VelTop; call Stream_and_Collide; call calc_dfc; 
      if ( (real(dfc)+aimag(dfc).lt.target_acc_Hz) .and. (errorprev.gt.target_acc_Hz ) ) then 
        nsteps(iovt) = it; errorprev = real(dfc)+aimag(dfc)
      endif  
      if(mod(it,int(nsteps_mx/10)) == 1) print *,'Liq',it,nsteps_mx,real(dfc),aimag(dfc)
    enddo
    dfc_ref(iovt) = dfc;
    B_Save = .true. 
    call compute_VelTop; call Stream_and_Collide; 
    open(unit = 20,file='ux_vs_y.asc1', action='write',status='unknown', position='append'); 
    do y = 1,ny 
      write (20,*) y*gridreso_SI*1e9,real(ux(Poij(1,y,1))), aimag(ux(Poij(1,y,1))); 
    enddo  
    close(unit = 20)
  enddo
  call DeAllocateAll

  rad_Sph = rad_Sph_SI * delta_f0 / delta_f0_SI; 
  nx = int(2. * rad_Sph * sizefac); 
  if (mod(nx,2).eq.0) nx = nx + 1;
  nz = nx; x_Sph = nx/2.; z_Sph = nz/2.; y_Sph = 0.
!  rad_Sph = 1.5 !!!!!!!!!!!!!!!!!!!!
  call system_clock(beginning,rate)
  call Count_and_Allocate
  call make_Lists
  do iovt = 1,novt
    B_Save = (i_calc_type.eq.j_elastic)
    call Set_om_dependent_variables(iovt); 
    print *,'nsteps(iovt)',nsteps(iovt)
    do it = 1,nsteps(iovt)
      call compute_VelTop; 
      if(i_calc_type.eq.j_BounceB) call sphere_dynamics
      call Stream_and_Collide; 
      call calc_dfc; 
      if(mod(it,int(nsteps(iovt)/10)) == 1) then 
        print *,it,nsteps(iovt),real(dfc- dfc_ref(iovt)),aimag(dfc- dfc_ref(iovt))
      endif  
      if (mod(it,int(nsteps(iovt)/100)) == 1) then 
        open(unit = 20,file='freq_ring_in.asc1',action='write',status='unknown', position='append'); 
        write (20,*)   ovtorder,it,nsteps(iovt),real(dfc- dfc_ref(iovt)),aimag(dfc- dfc_ref(iovt)); 
        close(unit = 20);   
      endif  
    enddo
    dfc_fin(iovt) = dfc - dfc_ref(iovt)
    B_Save = .true. 
    call compute_VelTop; call Stream_and_Collide; 
  enddo
  call system_clock(ending); print *,"Elapsed time is in seconds --->",real(ending - beginning)/real(rate)

  open(unit = 20,file='fin_freq.asc1',    action='write',status='unknown', position='append'); 
  write (20,*) cell_height_SI,ovtorder,real(dfc_fin)/ovtorder,aimag(dfc_fin)/ovtorder; close(unit = 20)
  
  if (i_calc_type .eq. j_BounceB) then
    do iPo = 1,nPoSph
      x = PoSphx(iPo)
      y = PoSphy(iPo)
      z = PoSphz(iPo)
      ux(Poij(x,y,z)) = 1
    enddo  
  endif  
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
  call DeAllocateAll
  deallocate (nsteps)
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
  omega = 2.*nu/delta**2.; 
  omega2by2 = omega**2/2; Zliq = -conjg(sqrt(rho0*complex(0,1)*omega*rho0*nu))
  tauinv_Liq_m_omega2by2 = tauinv_Liq - omega2by2
  tauinv_Sph_m_omega2by2 = tauinv_Sph - omega2by2
  omtauinv_Sph_p_Iomega = omtauinv_Sph + complex(0.,1.)*omega
  Iomega = complex(0.,1.)*omega
  K_Sph = nu * rho0 * omega * 1e4;
  nsteps_mx = int(ncycles / omega * (ny / delta)**2)
end subroutine Set_om_dependent_variables

subroutine Count_and_Allocate
  use declare; implicit none
  integer :: x,y,z
  logical :: in_Sph,in_Liq,on_WalLiq,on_Bor,on_3PL,on_Rim
  nPo = nx*ny*nz; nLi = nPo * (ndir + 1); 
  nPoLiq = 0; nPoSph = 0; nPoTop = 0; nPoWalLiq = 0; nPoWalSph = 0; nPoBor = 0; nPo3PL = 0; nPoRim = 0;  
  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        if ( (in_Liq(x,y,z)) .and. (y.gt.1) .and. (y.lt.ny) )   nPoLiq = nPoLiq + 1        
        if ( on_WalLiq(x,y,z) )                                 nPoWalLiq = nPoWalLiq + 1
        if ( (y.eq.1) .and. in_Sph(x,y,z) )                     nPoWalSph = nPoWalSph + 1
        if (in_Sph(x,y,z))      nPoSph = nPoSph + 1
        if (y.eq.ny)            nPoTop = nPoTop + 1
        if (on_Bor(x,y,z))      nPoBor = nPoBor + 1
        if (on_3PL(x,y,z))      nPo3PL = nPo3PL + 1
        if (on_Rim(x,y,z))      nPoRim = nPoRim + 1
      enddo
    enddo
  enddo
  allocate( h1(nLi),h2(nLi) )
  allocate( dr(nPo),ux(nPo),uy(nPo),uz(nPo) )
  if (i_calc_type.eq.j_elastic) allocate( dx(nPo),dy(nPo),dz(nPo) )
  allocate( uxTop(nPoTop), uzTop(nPoTop) )
  allocate( Rimux(nPoRim),Rimuy(nPoRim),Rimuz(nPoRim) )
  allocate( FxBor(nPoBor,0:ndir),FyBor(nPoBor,0:ndir),FzBor(nPoBor,0:ndir) )
  allocate( Fx3PLSph(nPo3PL,0:ndir),Fy3PLSph(nPo3PL,0:ndir),Fz3PLSph(nPo3PL,0:ndir) )
  allocate( Fx3PLWal(nPo3PL,0:ndir),Fy3PLWal(nPo3PL,0:ndir),Fz3PLWal(nPo3PL,0:ndir) )
  allocate( FxWalLiq(nPoWalLiq,1:5),FxWalSph(nPoWalSph,1:5) )
  allocate( dfc_ref(novt),dfc_fin(novt) )
  allocate( Poij(nx,ny,nz) )
  allocate( Rimj(nx,ny,nz) )
  allocate( Linj(nx,ny,nz,0:ndir) )
  allocate( Topj(nx,nz) )
  allocate( nLiBor(nPoBor),nLi3PLSph(nPo3PL),nLi3PLWal(nPo3PL) )
  allocate( PoLiqx(nPoLiq),PoLiqy(nPoLiq),PoLiqz(nPoLiq) )
  allocate( PoSphx(nPoSph),PoSphy(nPoSph),PoSphz(nPoSph) )
  allocate( PoTopx(nPoTop),PoTopy(nPoTop),PoTopz(nPoTop) )
  allocate( PoBorx(nPoBor),PoBory(nPoBor),PoBorz(nPoBor) )
  allocate( Po3PLx(nPo3PL),Po3PLy(nPo3PL),Po3PLz(nPo3PL) )
  allocate( PoRimx(nPoRim),PoRimy(nPoRim),PoRimz(nPoRim) )
  allocate( PoWalLiqx(nPoWalLiq),PoWalLiqy(nPoWalLiq),PoWalLiqz(nPoWalLiq) )
  allocate( PoWalSphx(nPoWalSph),PoWalSphy(nPoWalSph),PoWalSphz(nPoWalSph) )
  allocate( permu_Bor(ndir+1,nPoBor),permu_3PL(ndir+1,nPo3PL))
end subroutine Count_and_Allocate

subroutine DeAllocateAll
  use declare; implicit none
  deallocate(h1,h2); 
  deallocate(dr,ux,uy,uz);
  if (i_calc_type.eq.j_elastic) deallocate(dx,dy,dz)
  deallocate( uxTop, uzTop )
  deallocate( Rimux,Rimuy,Rimuz )
  deallocate( FxBor,FyBor,FzBor)
  deallocate( Fx3PLSph,Fy3PLSph,Fz3PLSph )
  deallocate( Fx3PLWal,Fy3PLWal,Fz3PLWal )
  deallocate( FxWalLiq,FxWalSph )
  deallocate( dfc_ref,dfc_fin )
  deallocate( Poij,Rimj,Linj,Topj )
  deallocate( nLiBor,nLi3PLSph,nLi3PLWal )
  deallocate( PoLiqx,PoLiqy,PoLiqz )
  deallocate( PoSphx,PoSphy,PoSphz )
  deallocate( PoTopx,PoTopy,PoTopz )
  deallocate( PoBorx,PoBory,PoBorz )
  deallocate( Po3PLx,Po3PLy,Po3PLz )
  deallocate( PoRimx,PoRimy,PoRimz )
  deallocate( PoWalLiqx,PoWalLiqy,PoWalLiqz )
  deallocate( PoWalSphx,PoWalSphy,PoWalSphz )
  deallocate( permu_Bor,permu_3PL)
end subroutine DeAllocateAll

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
      i = permu_Wal_BBa(iloc)
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
        i = permu_Wal_BBa(iloc)
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
      if (i_calc_type.eq.j_BounceB) then
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
      if (i_calc_type.eq.j_BounceB) then
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
      if (i_calc_type.eq.j_BounceB) then
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
      i = permu_Top_BBa(iloc)
      ibar = no_slip(i)
      cidotu = dirx(ibar)*uxTop(Topj(x,z)) + dirz(ibar)*uzTop(Topj(x,z))
      ht(i) = h1(Linj(x,y,z,ibar)) - 6 * wi(ibar) * cidotu
    enddo  
    do iloc = 1,14
      i = permu_Top_Liq(iloc)
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


