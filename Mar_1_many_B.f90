module declare
  implicit none
  integer, parameter :: ndir = 18
  integer, parameter :: dirx(0:18) = (/0,+1,-1, 0, 0, 0, 0,+1,-1,+1,-1, 0, 0,+1,-1,+1,-1, 0, 0/)
  integer, parameter :: dirz(0:18) = (/0, 0, 0, 0, 0,+1,-1, 0, 0,+1,-1,+1,-1, 0, 0,-1,+1,-1,+1/)
  integer, parameter :: diry(0:18) = (/0, 0, 0,+1,-1, 0, 0,+1,-1, 0, 0,+1,-1,-1,+1, 0, 0,+1,-1/)
  integer, parameter :: no_slip(0:18) = (/0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17/)
  real*8,  parameter :: wi(0:18) = (/ 1./3.,1./18,1./18.,1./18.,1./18.,1./18.,1./18., &
                           1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36./)
  real*8,parameter :: pi = 4.D0*DATAN(1.D0), rho0 = 1., nu = 1./6., tauinv_Liq = 2./(6.*nu + 1.) 
  real*8,parameter :: omtauinv_Liq = 1. - tauinv_Liq

  real*8,parameter :: Zq_SI = 8.8e6, nu_SI = 1e-6, rho0_SI = 1e3, f0_SI = 5e6, df_analy = 716.799 
  real*8,parameter :: delta_f0_SI = sqrt(2 * nu_SI/(2*pi*f0_SI));
  real*8           :: f_SI, cell_height_SI, grid_reso_SI, rad_Sph_SI, rho_Sph_SI
  real*8           :: truncation_by_diameter, cell_diam_by_sphere_diam
  real*8           :: delta_f0, f0, rad_Sph, rho_Sph
  real*8           :: ncycles_mx, delta, omega, omega2by2 
  real*8 :: mass_Sph,target_acc_Hz
  complex*16 :: Iomega, tauinv_Liq_m_omega2by2, FxTotOnWal, dfc, ZL, Zq, Zliq
  complex*16 :: dfc_ref,dfc_fin,dfc_upb,u_upb,u_upb_ref 

  integer :: nx,ny,nz,nLi,nPo,nPoLiq,nPoSph,nPoTop,nPoWalLiq,nPoBor,nPo3PL
  integer :: it,nsteps,nsteps_mx,error,i_BBa_simple_inpol,ovtorder
  integer :: beginning,rate,ending
  integer, parameter :: n_Sph_max = 10
  integer :: i_Sph,n_Sph
  real*8 :: x_Sph(n_Sph_max), y_Sph(n_Sph_max), z_Sph(n_Sph_max),xfrac_Sph(n_Sph_max),zfrac_Sph(n_Sph_max) 
!  integer :: x_Sph,y_Sph,z_Sph

  real*8    , allocatable ::  xBor(:,:), yBor(:,:), zBor(:,:), qBor(:,:)
  complex*16, allocatable :: dr(:),ux(:),uy(:),uz(:),h1(:),h2(:)
  complex*16, allocatable :: uxTop(:), uzTop(:)
  complex*16, allocatable :: FxBor(:,:),FyBor(:,:),FzBor(:,:)
  complex*16, allocatable :: uxBor(:,:),uyBor(:,:),uzBor(:,:)
  complex*16, allocatable :: Fx3PLOnSph(:,:), Fy3PLOnSph(:,:),Fz3PLOnSph(:,:)
  complex*16, allocatable :: Fx3PLOnWal(:,:), FxLiqOnWal(:,:)
  complex*16, allocatable :: Fx(:,:,:),Fy(:,:,:),Fz(:,:,:)
    
  integer, allocatable :: Poij(:,:,:),Linj(:,:,:,:),Topj(:,:)
  integer, allocatable :: nLiBor(:),nLi3PLSph(:),nLi3PLWal(:)
  integer, allocatable :: PoLiqx(:),PoLiqy(:),PoLiqz(:)
  integer, allocatable :: PoSphx(:),PoSphy(:),PoSphz(:)
  integer, allocatable :: PoTopx(:),PoTopy(:),PoTopz(:)
  integer, allocatable :: PoBorx(:),PoBory(:),PoBorz(:)
  integer, allocatable :: Po3PLx(:),Po3PLy(:),Po3PLz(:)
  integer, allocatable :: PoWalLiqx(:),PoWalLiqy(:),PoWalLiqz(:)
  integer, allocatable :: permu_Bor(:,:),permu_3PL(:,:)
  integer, parameter :: permu_Top_BBa(1:5 ) = (/        4,      8,        13,12,            18/) !page87
  integer, parameter :: permu_Top_Liq(1:14) = (/0,1,2,3,  5,6,7,  9,10,11,      14,15,16,17   /)
  integer, parameter :: permu_Wal_BBA(1:5 ) = (/      3,      7,       14,      11,      17   /)
  integer, parameter :: permu_Wal_Liq(1:14) = (/0,1,2,  4,5,6,  8,9,10,   12,13,   15,16,   18/)
end module declare

program FT_LBM
  use declare
  implicit none
  complex * 16 :: dfc_prev 
  logical :: B_accuracy_reached
  integer iovt
  
  open(unit = 20,file='freq_ring_in.asc1',    action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='freq_ring_in_upb.asc1',action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='ux_vs_y.asc1',         action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='uxuy_3D.asc1',         action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='sim_time_steps.asc1',  action='write',status='unknown'); close(unit = 20,status='delete');   

  open(unit = 21,file='sim_config_B.asc1',      action='read',status='unknown'); 
  read(21,*) target_acc_Hz
  read(21,*) rad_Sph_SI
  read(21,*) cell_height_SI
  read(21,*) grid_reso_SI
  print *,'grid_reso_SI',grid_reso_SI
  read(21,*) cell_diam_by_sphere_diam
  print *,'cell_diam_by_sphere_diam',cell_diam_by_sphere_diam
  read(21,*) truncation_by_diameter
  print *,'truncation_by_diameter',truncation_by_diameter
  read(21,*) iovt
  print *,'iovt',iovt
  read(21,*) n_Sph
  do i_Sph = 1,n_Sph 
    read(21,*) xfrac_Sph(i_Sph),zfrac_Sph(i_Sph)
  enddo
 
  print *,'iovt',   iovt
  close(unit = 21);   
 

!  target_acc_Hz = 0.1
!  rad_Sph_SI = 5e-9
!  cell_height_SI = 20e-9
!  grid_reso_SI  = 1.5e-9
!  cell_diam_by_sphere_diam = 1.5
!  truncation_by_diameter = 0.5
!  iovt = 1
  i_BBa_simple_inpol = 1
  

  ncycles_mx = 20; 
  ny = int(cell_height_SI/grid_reso_SI)
  delta_f0 = delta_f0_SI/grid_reso_SI
  f0 = 2.*nu/delta_f0**2 / (2. * pi)

  rad_Sph = 0; nx = 1; nz = 1; nPoTop = nx * nz; 

  call Count_and_Allocate
  call make_Lists
  call Set_om_dependent_variables(iovt); 
  nsteps = nsteps_mx
  call init_u0; 
  do it = 1,nsteps_mx
    call compute_VelTop; 
    call Stream_and_Collide; 
    call calc_dfc; 
    if(mod(it,int(nsteps/5)) == 1) then 
      print *,'Ref State',it,nsteps,real(dfc),aimag(dfc)
    endif  
  enddo
  dfc_ref = dfc;
  call calc_u_upb(u_upb_ref)
  call DeAllocateAll
  rad_Sph = rad_Sph_SI * delta_f0 / delta_f0_SI; 
  nx = int(2. * rad_Sph * cell_diam_by_sphere_diam); 
  if (mod(nx,2).eq.0) nx = nx + 1;
  nz = nx;
  do i_Sph = 1,n_Sph 
    x_Sph(i_Sph) = xfrac_Sph(i_Sph)*nx; 
    z_Sph(i_Sph) = zfrac_Sph(i_Sph)*nz; 
    y_Sph(i_Sph) = rad_Sph * (1 - 2. * truncation_by_diameter) + 1
  enddo  
  print *, nx,ny,nz,x_Sph(1),y_Sph(1),z_Sph(1) 
!  x_Sph = int(nx/2.)+1; z_Sph = int(nz/2.)+1; y_Sph = rad_Sph * (1 - 2. * truncation_by_diameter) + 1

  call system_clock(beginning,rate)
  call Count_and_Allocate
  call make_Lists
  print *,'Lists done'
  call Set_om_dependent_variables(iovt); 
  call init_u0; 
  it = 0  
  B_accuracy_reached = .false.
  do while (it < nsteps_mx .and. .not.B_accuracy_reached )
    it = it + 1  
    call compute_VelTop; 
    call Stream_and_Collide; 
    call calc_dfc;
    if ( abs(dfc_prev - dfc) < target_acc_Hz ) then 
      print *,'limit reached',it
      B_accuracy_reached = .true.
    endif  
    dfc_prev = dfc
    if(mod(it,10) == 1) then 
      print "(I5,I7,F16.5,F16.5,F12.5,F12.5)",it,nsteps,real(dfc),aimag(dfc),real(dfc_upb),aimag(dfc_upb)
      open(unit = 20,file='freq_ring_in.asc1',action='write',status='unknown', position='append'); 
      write (20,*)   it,real(dfc-dfc_ref),aimag(dfc-dfc_ref); 
      close(unit = 20);   
      open(unit = 20,file='freq_ring_in_upb.asc1',action='write',status='unknown', position='append'); 
      write (20,*) it,real(dfc_upb),aimag(dfc_upb); 
      close(unit = 20);   
      call write_uxuy
    endif  
  enddo
  dfc_fin = dfc - dfc_ref
  call compute_VelTop; call Stream_and_Collide; 
  call system_clock(ending); print *,"Elapsed time [secs] : ",real(ending - beginning)/real(rate)
  open(unit = 20,file='sim_time_steps.asc1',action='write'); 
  write(20,"(I5,I5,I8,I6)") nx,ny,nPo,nsteps
  write(20,"(F16.8)") real(ending - beginning)/real(rate)
  write(20,"(E17.8,E17.8,E17.8)") real(dfc_fin), aimag(dfc_fin), grid_reso_SI
  close(unit = 20);   
  call Plot_Fx
  call DeAllocateAll
end program FT_LBM

subroutine init_u0
  use declare; implicit none
  integer :: x,y,z,iPo,iLi
  do iPo = 1,nPo
    dr(iPo) = 0
    ux(iPo) = 0
    uy(iPo) = 0
    uz(iPo) = 0
  enddo  
  do iPo = 1,nPoSph
    x = PoSphx(iPo)
    y = PoSphy(iPo)
    z = PoSphz(iPo)
    ux(Poij(x,y,z)) = 1
  enddo  
  do iLi = 1,nLi
    h1(iLi) = 0
  enddo
end subroutine init_u0

logical function in_Sph(x,y,z) 
  use declare; implicit none
  integer :: x,y,z
  in_Sph = .false.
  do i_Sph = 1,n_Sph
    if ( (y.ge.1) .and. ((x-x_Sph(i_Sph))**2 + (y-y_Sph(i_Sph))**2 + (z-z_Sph(i_Sph))**2 .lt. rad_Sph**2) ) then 
!    if ( (y.ge.1) .and. ((x-x_Sph)**2 + (y-y_Sph)**2 + (z-z_Sph)**2 .lt. rad_Sph**2) ) then 
      in_Sph = .true.
    endif  
  enddo
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

subroutine Set_om_dependent_variables(iovt)
  use declare; implicit none
  integer :: iovt
  ovtorder = 2 * iovt - 1; f_SI = ovtorder* f0_SI; 
  delta = delta_f0 / sqrt(float(ovtorder)); 
  omega = 2.*nu/delta**2.; 
  omega2by2 = omega**2/2; Zliq = sqrt(rho0*complex(0,1)*omega*rho0*nu)
  tauinv_Liq_m_omega2by2 = tauinv_Liq - omega2by2
  Iomega = complex(0.,1.)*omega
  nsteps_mx = int(ncycles_mx / omega * (ny / delta)**2)
end subroutine Set_om_dependent_variables

subroutine Count_and_Allocate
  use declare; implicit none
  integer :: x,y,z
  logical :: in_Sph,in_Liq,on_WalLiq,on_Bor,on_3PL
  nPo = nx*ny*nz; nLi = nPo * (ndir + 1); 
  nPoLiq = 0; nPoSph = 0; nPoBor = 0; nPoTop = 0; nPoWalLiq = 0; nPo3PL = 0; 
  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        if ( (in_Liq(x,y,z)) .and. (y.gt.1) .and. (y.lt.ny) ) &
                                nPoLiq = nPoLiq + 1        
        if ( on_WalLiq(x,y,z) ) nPoWalLiq = nPoWalLiq + 1
        if (in_Sph(x,y,z))      nPoSph = nPoSph + 1
        if (y.eq.ny)            nPoTop = nPoTop + 1
        if (on_Bor(x,y,z))      nPoBor = nPoBor + 1
        if (on_3PL(x,y,z))      nPo3PL = nPo3PL + 1
      enddo
    enddo
  enddo
  allocate( h1(nLi),h2(nLi) )
  allocate( dr(nPo),ux(nPo),uy(nPo),uz(nPo) )
  allocate( uxTop(nPoTop), uzTop(nPoTop) )
  allocate( FxBor(nPoBor,0:ndir),FyBor(nPoBor,0:ndir),FzBor(nPoBor,0:ndir))
  allocate( qBor(nPoBor,0:ndir), xBor(nPoBor,0:ndir),yBor(nPoBor,0:ndir),zBor(nPoBor,0:ndir) )
  allocate( uxBor(nPoBor,0:ndir),uyBor(nPoBor,0:ndir),uzBor(nPoBor,0:ndir) )
  allocate( Fx3PLOnSph(nPo3PL,0:ndir), Fy3PLOnSph(nPo3PL,0:ndir), Fz3PLOnSph(nPo3PL,0:ndir) )
  allocate( Fx3PLOnWal(nPo3PL,0:ndir), FxLiqOnWal(nPoWalLiq,1:5) )
  allocate( Poij(nx,ny,nz),Topj(nx,nz),Linj(nx,ny,nz,0:ndir) )
  allocate( nLiBor(nPoBor),nLi3PLSph(nPo3PL),nLi3PLWal(nPo3PL) )
  allocate( PoLiqx(nPoLiq),PoLiqy(nPoLiq),PoLiqz(nPoLiq) )
  allocate( PoSphx(nPoSph),PoSphy(nPoSph),PoSphz(nPoSph) )
  allocate( PoTopx(nPoTop),PoTopy(nPoTop),PoTopz(nPoTop) )
  allocate( PoBorx(nPoBor),PoBory(nPoBor),PoBorz(nPoBor) )
  allocate( Po3PLx(nPo3PL),Po3PLy(nPo3PL),Po3PLz(nPo3PL) )
  allocate( PoWalLiqx(nPoWalLiq),PoWalLiqy(nPoWalLiq),PoWalLiqz(nPoWalLiq) )
  allocate( permu_Bor(nPoBor,ndir+1),permu_3PL(nPo3PL,ndir+1))
  allocate( Fx(nx,ny,nz),Fy(nx,ny,nz),Fz(nx,ny,nz) )
end subroutine Count_and_Allocate

subroutine DeAllocateAll
  use declare; implicit none
  deallocate(h1,h2); 
  deallocate(dr,ux,uy,uz);
  deallocate( uxTop, uzTop )
  deallocate( FxBor,FyBor,FzBor)
  deallocate( qBor,xBor,yBor,zBor )  
  deallocate( uxBor,uyBor,uzBor )  
  deallocate( Fx3PLOnSph,Fy3PLOnSph,Fz3PLOnSph )
  deallocate( Fx3PLOnWal,FxLiqOnWal )
  deallocate( Poij,Linj,Topj )
  deallocate( nLiBor,nLi3PLSph,nLi3PLWal )
  deallocate( PoLiqx,PoLiqy,PoLiqz )
  deallocate( PoSphx,PoSphy,PoSphz )
  deallocate( PoTopx,PoTopy,PoTopz )
  deallocate( PoBorx,PoBory,PoBorz )
  deallocate( Po3PLx,Po3PLy,Po3PLz )
  deallocate( PoWalLiqx,PoWalLiqy,PoWalLiqz )
  deallocate( permu_Bor,permu_3PL)
  deallocate( Fx,Fy,Fz ); 
end subroutine DeAllocateAll

subroutine make_Lists
  use declare; implicit none
  integer :: x,y,z,xmd,ymd,zmd,i,i_next
  logical :: in_Sph,in_Liq,on_WalLiq,on_Bor,on_3PL
  nPo = 0; nLi = 0; nPoLiq = 0; nPoSph = 0; nPoTop = 0; nPoWalLiq = 0; 
  nPoBor = 0; nPo3PL = 0; 
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
              permu_Bor(nPoBor,nLiBor(nPoBor)) = i
              call find_q(x,y,z,xmd,ymd,zmd,qBor(nPoBor,i),xBor(nPoBor,i),yBor(nPoBor,i),zBor(nPoBor,i))
            endif
          enddo
          i_next = nLiBor(nPoBor)
          do i = 0,ndir
            xmd = mod(nx+x-dirx(i)-1,nx)+1
            ymd = y-diry(i)
            zmd = mod(nz+z-dirz(i)-1,nz)+1
            if ((ymd.ge.1) .and. (.not.in_Sph(xmd,ymd,zmd)) ) then 
              i_next = i_next + 1
              permu_Bor(nPoBor,i_next) = i
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
              permu_3PL(nPo3PL,nLi3PLSPh(nPo3PL)) = i
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
              permu_3PL(nPo3PL,i_next) = i
            endif
          enddo
          do i = 0,ndir
            xmd = mod(nx+x-dirx(i)-1,nx)+1
            ymd = y-diry(i)
            zmd = mod(nz+z-dirz(i)-1,nz)+1
            if ( ( .not. in_Sph(xmd,ymd,zmd)) .and. (ymd.ge.1) ) then 
              i_next = i_next + 1
              permu_3PL(nPo3PL,i_next) = i
            endif
          enddo
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
  h2(Linj(x,y,z,i)) = tauinv_Liq_m_omega2by2 * heq - Iomega * ht(i)   !omtauinv = 0
end subroutine relax_Liq

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
  dr(Poij(x,y,z)) = drloc
  ux(Poij(x,y,z)) = uxloc
  uy(Poij(x,y,z)) = uyloc
  uz(Poij(x,y,z)) = uzloc
end subroutine calc_dr_uxyz


subroutine Stream_and_Collide
  use omp_lib
  use declare; implicit none
  integer :: x,y,z,i,iloc,ibar,xmd,ymd,zmd,xpd,ypd,zpd,iPo,iLi
  complex*16 :: ht(0:ndir),drloc,uxloc,uyloc,uzloc,cidotu
  real* 8 :: q
  logical in_Liq

!  !$omp parallel private(iPo,i,x,y,z,drloc,uxloc,uyloc,xmd,ymd,zmd,uzloc,ht)
!    !$omp do
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
!    !$omp end do
!  !$omp end parallel

  do iPo = 1,nPoWalLiq
    x = PoWalLiqx(iPo)
    y = PoWalLiqy(iPo)
    z = PoWalLiqz(iPo)
    do iloc = 1,5
      ibar = permu_Wal_BBa(iloc)
      i = no_slip(ibar)
      ht(ibar) = h1(Linj(x,y,z,i)) - 6 * wi(i) * dirx(i)  ! page 180 
      FxLiqOnWal(iPo,iloc) = (ht(ibar) + h1(Linj(x,y,z,i))) * dirx(i)   ! page 217
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

  do iPo = 1,nPoBor
    x = PoBorx(iPo)
    y = PoBory(iPo)
    z = PoBorz(iPo)
    do iloc = 1,nLiBor(iPo)
      ibar = permu_Bor(iPo,iloc)
      i = no_slip(ibar)
      cidotu = dirx(i)
      if (i_BBa_simple_inpol == 1) then
        ht(ibar) = h1(Linj(x,y,z,i)) - 6 * wi(i) * cidotu
      endif
      if (i_BBa_simple_inpol == 2) then
        xpd = mod(nx+x+dirx(i)-1,nx)+1
        ypd = y+diry(i)
        zpd = mod(nz+z+dirz(i)-1,nz)+1
        q = qBor(iPo,i)
        if (q < 1./2.) then
          if ( (in_Liq(xpd,ypd,zpd)) ) then 
            ht(ibar) = 2*q * h1(Linj(x,y,z,i)) + (1-2*q) * h1(Linj(xpd,ypd,zpd,i)) - 6 * wi(i) * cidotu
          else
            ht(ibar) = h1(Linj(x,y,z,i)) - 6 * wi(i) * cidotu
          endif
        endif
        if (q >= 1./2.) then 
          ht(ibar) = 1./(2*q) * h1(Linj(x,y,z,i)) + (2*q-1)/(2*q) * h1(Linj(x,y,z,i)) - 6 * wi(i) * cidotu
        endif
      endif
      FxBor(iPo,iloc) = (ht(ibar) + h1(Linj(x,y,z,i))) * dirx(i)
      FyBor(iPo,iloc) = (ht(ibar) + h1(Linj(x,y,z,i))) * diry(i)
      FzBor(iPo,iloc) = (ht(ibar) + h1(Linj(x,y,z,i))) * dirz(i)
    enddo  
    do iloc = nLiBor(iPo)+1,ndir+1
      i = permu_Bor(iPo,iloc)
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
      ibar = permu_3PL(iPo,iloc)
      i = no_slip(ibar)
      cidotu = dirx(i)
      if (i_BBa_simple_inpol == 1) then
        ht(ibar) = h1(Linj(x,y,z,i)) - 6 * wi(i) * cidotu
      endif
      if (i_BBa_simple_inpol == 2) then
        xpd = mod(nx+x+dirx(i)-1,nx)+1
        ypd = y+diry(i)
        zpd = mod(nz+z+dirz(i)-1,nz)+1
        q = qBor(iPo,i)
        if (q < 1./2.) then
          if ( (in_Liq(xpd,ypd,zpd)) ) then 
            ht(ibar) = 2*q * h1(Linj(x,y,z,i)) + (1-2*q) * h1(Linj(xpd,ypd,zpd,i)) - 6 * wi(i) * cidotu
          else
            ht(ibar) = h1(Linj(x,y,z,i)) - 6 * wi(i) * cidotu
          endif
        endif
        if (q >= 1./2.) then 
          ht(ibar) = 1./(2*q) * h1(Linj(x,y,z,i)) + (2*q-1)/(2*q) * h1(Linj(x,y,z,i)) - 6 * wi(i) * cidotu
        endif
      endif
      Fx3PLOnSph(iPo,iloc) = (ht(ibar) + h1(Linj(x,y,z,i))) * dirx(i)
      Fy3PLOnSph(iPo,iloc) = (ht(ibar) + h1(Linj(x,y,z,i))) * diry(i)
      Fz3PLOnSph(iPo,iloc) = (ht(ibar) + h1(Linj(x,y,z,i))) * dirz(i)
    enddo  
    do iloc = nLi3PLSPh(iPo) + 1, nLi3PLSPh(iPo) + nLi3PLWal(iPo)
      ibar = permu_3PL(iPo,iloc)
      i = no_slip(ibar)
      cidotu = dirx(i)
      if (i_BBa_simple_inpol == 1) then
        ht(ibar) = h1(Linj(x,y,z,i)) - 6 * wi(i) * cidotu
      endif
      if (i_BBa_simple_inpol == 2) then
        xpd = mod(nx+x+dirx(i)-1,nx)+1
        ypd = y+diry(i)
        zpd = mod(nz+z+dirz(i)-1,nz)+1
        q = qBor(iPo,i)
        if (q < 1./2.) then
          if ( (in_Liq(xpd,ypd,zpd)) ) then 
            ht(ibar) = 2*q * h1(Linj(x,y,z,i)) + (1-2*q) * h1(Linj(xpd,ypd,zpd,i)) - 6 * wi(i) * cidotu
          else
            ht(ibar) = h1(Linj(x,y,z,i)) - 6 * wi(i) * cidotu
          endif
        endif
        if (q >= 1./2.) then 
          ht(ibar) = 1./(2*q) * h1(Linj(x,y,z,i)) + (2*q-1)/(2*q) * h1(Linj(x,y,z,i)) - 6 * wi(i) * cidotu
        endif
      endif
      Fx3PLonWal(iPo,iloc) = (ht(ibar) + h1(Linj(x,y,z,i))) * dirx(i)
    enddo  
    do iloc = nLi3PLSPh(iPo) + nLi3PLWal(iPo) + 1,ndir+1
      i = permu_3PL(iPo,iloc)
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

  do iLi = 1,nLi
    h1(iLi) = h2(iLi)
  enddo
end subroutine Stream_and_Collide


subroutine calc_dfc
  use declare; implicit none
  integer :: iPo,iloc
  complex*16 :: FxTotOnWalLiq,FxTotOnWalSph,dfcliq,dfcSph,dfcMas
  FxTotOnWal = 0.
  FxTotOnWalLiq = 0.
  FxTotOnWalSph = 0.
  
  do iPo = 1,nPoWalLiq
    do iloc = 1,5
      FxTotOnWal    = FxTotOnWal    + FxLiqOnWal(iPo,iloc)     
      FxTotOnWalLiq = FxTotOnWalLiq + FxLiqOnWal(iPo,iloc)     
    enddo
  enddo
  do iPo = 1,nPo3PL
    do iloc = 1, nLi3PLSph(iPo)
      FxTotOnWal    = FxTotOnWal    + Fx3PLOnSph(iPo,iloc)
      FxTotOnWalSph = FxTotOnWalSph + Fx3PLOnSph(iPo,iloc)
    enddo
    do iloc = nLi3PLSph(iPo) + 1, nLi3PLSph(iPo) + nLi3PLWal(iPo) 
      FxTotOnWal    = FxTotOnWal    + Fx3PLOnWal(iPo,iloc)
      FxTotOnWalLiq = FxTotOnWalLiq + Fx3PLOnWal(iPo,iloc)
    enddo
  enddo
  do iPo = 1,nPoBor
    do iloc = 1,nLiBor(iPo)
      FxTotOnWal    = FxTotOnWal    + FxBor(iPo,iloc)
      FxTotOnWalSph = FxTotOnWalSph + FxBor(iPo,iloc)
    enddo
  enddo

  mass_Sph = rho0 * rad_Sph**3. * 4.*pi/3 * &
      (truncation_by_diameter-1.)**2 * (1 + 2*truncation_by_diameter)  !4/3 (-1 + al)^2 (1 + 2 al) \[Pi]

  ZL = -(FxTotOnWal)/(nx*nz) - Iomega * mass_Sph /(nx*nz)
  Zq = (rho0/rho0_SI)*(nu/nu_SI)*(delta_f0_SI/delta_f0)*Zq_SI
  dfc    = - f0_SI * complex(0,1)/(pi*Zq) * ZL
  dfcLiq = - f0_SI * complex(0,1)/(pi*Zq) * (-(FxTotOnWalLiq)/(nx*nz))
  dfcSph = - f0_SI * complex(0,1)/(pi*Zq) * (-(FxTotOnWalSph)/(nx*nz))
  dfcMas = - f0_SI * complex(0,1)/(pi*Zq) * (- Iomega * mass_Sph /(nx*nz))
  call calc_u_upb(u_upb)
  ZL = Zliq * (u_upb / u_upb_ref - 1) - Iomega * mass_Sph /(nx*nz)
  dfc_upb    = - f0_SI * complex(0,1)/(pi*Zq) * ZL

!  if(mod(it,int(nsteps/10)) == 1) then 
!    print *,'dfc Liq',it,real(dfcLiq-dfc_ref),aimag(dfcLiq-dfc_ref)
!    print *,'dfc Sph',it,real(dfcSph-dfc_ref),aimag(dfcSph-dfc_ref)
!    print *,'dfc Liq + Sph',it,real(dfcLiq  + dfcSph-dfc_ref),aimag(dfcLiq  + dfcSph-dfc_ref)
!    print *,'dfc Mas',it,real(dfcMas),aimag(dfcMas)
!  endif 
end subroutine calc_dfc

subroutine find_q(x,y,z,xmd,ymd,zmd,q,xw,yw,zw)
  use declare; implicit none
  integer :: x,y,z,xmd,ymd,zmd,ibisec 
  real*8 :: xi,yi,zi,qprev,qnew,q,qc,xw,yw,zw
  qc = 1./2.
  qprev = 0
  do ibisec = 1,14
    xi = x + (xmd - x) * qc
    yi = y + (ymd - y) * qc
    zi = z + (zmd - z) * qc
    do i_Sph = 1,n_Sph 
      if ( (xi-x_Sph(i_Sph))**2 + (yi-y_Sph(i_Sph))**2 + (zi-z_Sph(i_Sph))**2 .lt. rad_Sph**2 ) then 
!      if ( (xi-x_Sph)**2 + (yi-y_Sph)**2 + (zi-z_Sph)**2 .lt. rad_Sph**2 ) then 
        qnew = qc - abs(qc - qprev) / 2.
      else
        qnew = qc + abs(qc - qprev) / 2.
      endif
    enddo  
    qprev = qc
    qc = qnew
    xw = xi; yw = yi; zw = zi; q = qnew
  enddo
end subroutine find_q

subroutine Plot_Fx
  use declare; implicit none
  integer :: x,y,z,i,iPo 
  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        Fx(nx,ny,nz) = 0; Fy(nx,ny,nz) = 0; Fz(nx,ny,nz) = 0;
      enddo
    enddo
  enddo
  do iPo = 1,nPoBor
    x = PoBorx(iPo)
    y = PoBory(iPo)
    z = PoBorz(iPo)
    do i = 0,ndir
      Fx(x,y,z) = Fx(x,y,z) + FxBor(iPo,i)
      Fy(x,y,z) = Fy(x,y,z) + FyBor(iPo,i)
      Fy(x,y,z) = Fy(x,y,z) + FyBor(iPo,i)
    enddo
  enddo  
  do iPo = 1,nPoWalLiq
    x = PoWalLiqx(iPo)
    y = 1
    z = PoWalLiqz(iPo)
    do i = 1,5
      Fx(x,1,z) = Fx(x,1,z) + FxLiqOnWal(iPo,i)
    enddo
!    print *,'Wal',iPo,x,y,z,real(Fx(x,1,z))
  enddo  

  do iPo = 1,nPo3PL
    x = Po3PLx(iPo)
    y = 1
    z = Po3PLz(iPo)
    do i = 0,ndir
      Fx(x,1,z) = Fx(x,1,z) + Fx3PLOnWal(iPo,i)
    enddo
!    print *,'3PL',iPo,x,y,z,real(Fx(x,1,z))
  enddo  

  open(unit = 20,file='FxFy_3D.asc1',action='write')
  write(20,"(I5,I5,I5)")nx,ny,nz
  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        write(20,"(I5,I5,I5,F12.8,F12.8,F12.8,F12.8,F12.8,F12.8)") x,y,z, &
                                   real(Fx(x,y,z)),aimag(Fx(x,y,z)),&
                                   real(Fy(x,y,z)),aimag(Fy(x,y,z)),&
                                   real(Fz(x,y,z)),aimag(Fz(x,y,z))
      enddo
    enddo
  enddo
  close(unit = 20)
end subroutine Plot_Fx

subroutine write_uxuy
  use declare; implicit none
  integer :: x,y
  open(unit = 20,file='uxuy_3D.asc1',action='write')
  write(20,"(I5,I5,I5)")nx,ny,nz
  do x = 1,nx
    do y = 1,ny
      write(20,"(I5,I5,F12.8,F12.8,F12.8,F12.8)") x,y, &
        real(ux(Poij(x,y,int(nz/2.)+1))),aimag(ux(Poij(x,y,int(nz/2.)+1))),&
        real(uy(Poij(x,y,int(nz/2.)+1))),aimag(uy(Poij(x,y,int(nz/2.)+1)))
    enddo
  enddo
  close(unit = 20)
end subroutine write_uxuy

subroutine calc_u_upb(u_upb_loc)
  use declare; implicit none
  integer :: x,z
  complex * 16 :: u_upb_loc
  u_upb_loc = 0
  do x = 1,nx
    do z = 1,nz
      u_upb_loc = u_upb_loc + uxTop(Topj(x,z)) 
    enddo  
  enddo  
  u_upb_loc = u_upb_loc / nPoTop
end subroutine calc_u_upb
