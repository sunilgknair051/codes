module declare
  implicit none
  real*8,parameter :: pi = 4.D0*DATAN(1.D0)
  complex*16,parameter :: iu = complex(0,1)

  integer, parameter :: ndir = 18
  integer, parameter :: dirx(0:18) = (/0,+1,-1, 0, 0, 0, 0,+1,-1,+1,-1, 0, 0,+1,-1,+1,-1, 0, 0/), &
                        dirz(0:18) = (/0, 0, 0, 0, 0,+1,-1, 0, 0,+1,-1,+1,-1, 0, 0,-1,+1,-1,+1/), &
                        diry(0:18) = (/0, 0, 0,+1,-1, 0, 0,+1,-1, 0, 0,+1,-1,-1,+1, 0, 0,+1,-1/), &
                        no_slip(0:18) = (/0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17/)
  real*8,  parameter :: wi(0:18) = (/ 1./3.,1./18,1./18.,1./18.,1./18.,1./18.,1./18., &
                           1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36./)
  real*8,parameter :: rho0 = 1., nu = 1./6., tauinv = 2./(6.*nu + 1.), eta0 = rho0 * nu, omtauinv = 1. - tauinv, &
                      Zq_SI = 8.8e6, nu_SI = 1e-6, rho0_SI = 1e3, f0_SI = 5e6, df_analy = 716.799, & 
                      delta_f0_SI = sqrt(2 * nu_SI/(2*pi*f0_SI));
  real*8           :: f_SI, cell_height_SI, grid_reso_SI, rad_Sph_SI, rho_Sph_SI, &
                      rad_Sph, rho_Sph, truncation_by_diameter, sphere_diam_by_cell_diam, &
                      delta_f0, f0, ncycles_mx, & 
                      delta, omega, omega2by2, tauinv_m_omega2by2 
  real*8 :: target_ratio, dif_ini

  complex*16 :: Iom,Zliq,uzTop_Avg,uxTop_Avg,uxTop_Avg_ref
  complex*16 :: dfc,dfc_ref,dfc_noMassLiq,dfc_noMassLiq_ref,dfc_stress,dfc_stress_ref,dfc_upb,dfc_upb_ref

  integer :: nx,ny,nz,nLi,nPo,nPoLiq,nPoSph,nPoTop,nPoWalLiq,nPoBor,nPo3PL, &
             it,nsteps,nsteps_mx,iovt,ovtorder, &
             beginning,rate,ending
  integer, parameter :: n_Sph_max = 10
  integer :: i_Sph,n_Sph
  real*8 :: x_Sph(n_Sph_max), y_Sph(n_Sph_max), z_Sph(n_Sph_max),xfrac_Sph(n_Sph_max),zfrac_Sph(n_Sph_max) 
  real*8 :: Mass_Liq, Qdot_Liq

  real*8    , allocatable :: y_arr_nm(:) 
  real*8    , allocatable :: xBor(:,:), yBor(:,:), zBor(:,:)
  complex*16, allocatable :: dr(:),ux(:),uy(:),uz(:),h1(:),h2(:)
  complex*16, allocatable :: uxTop(:), uzTop(:)
    
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
  logical :: B_save_u, B_dif_ini_found
  complex*16 :: FxLiqOnSphe, FxLiqOnWall,ZL
  real*8, parameter :: rho_Liq_SI = 1e3
  real*8 :: mass_Sph_SI, phase_diss

end module declare

program FT_LBM
  use declare
  implicit none
  complex * 16 :: dfc_prev,dfc_prev_prev 
  logical :: B_accuracy_reached
  B_save_u = .false.
  open(unit = 20,file='freq_ring_in.asc',    action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='freq_ring_in_stress.asc',    action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='freq_ring_in_upb.asc',    action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='freq_ring_in_noMassLiq.asc',    action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='freq_ring_in__noMassLiq.asc',    action='write',status='unknown'); close(unit = 20,status='delete');   
!  open(unit = 20,file='freq_ring_in_sqr_phase_diss.asc',action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='ux_vs_y.asc',         action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='uxuy_3D.asc',         action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='sim_time_steps.asc',  action='write',status='unknown'); close(unit = 20,status='delete');   

  open(unit = 21,file='sim_config.asc',      action='read',status='unknown'); 
  read(21,*) target_ratio
  read(21,*) rad_Sph_SI
  read(21,*) cell_height_SI
  read(21,*) grid_reso_SI
  read(21,*) sphere_diam_by_cell_diam
  read(21,*) truncation_by_diameter
  read(21,*) iovt
  read(21,*) n_Sph
  do i_Sph = 1,n_Sph 
    read(21,*) xfrac_Sph(i_Sph),zfrac_Sph(i_Sph)
  enddo
  close(unit = 21);   
 
  ncycles_mx = 10; 
  ny = int(cell_height_SI/grid_reso_SI)
  rad_Sph = 0; nx = 1; nz = 1; 
  call Count_and_Allocate
  call make_Lists
  call Set_om_dependent_variables; 
  nsteps = nsteps_mx
  call init_u0; 
  do it = 1,nsteps_mx
    call compute_VelTop; 
    call Stream_and_Collide; 
    if(mod(it,int(nsteps/5)) == 1) then 
      call calc_dfc; 
      print "(I5,F12.3,F12.3,F12.3,F12.3,F12.3,F12.3)",it,real(dfc)       ,aimag(dfc),&
                                                          real(dfc_stress),aimag(dfc_stress),&
                                                          real(dfc_upb),   aimag(dfc_upb)
    endif  
  enddo
  dfc_ref = dfc;
  dfc_stress_ref = dfc_stress 
  dfc_upb_ref = -716.8 + iu*716.8;
  dfc_noMassLiq_ref = dfc_noMassLiq;

  print *, dfc_upb_ref
  uxTop_Avg_ref = uxTop_Avg;
  call DeAllocateAll
  rad_Sph = rad_Sph_SI * delta_f0 / delta_f0_SI; 
  nx = int(2. * rad_Sph / sphere_diam_by_cell_diam); 
  print *, 'nx',nx
  print *,'rad_Sph',rad_Sph,sphere_diam_by_cell_diam
  if (mod(nx,2).eq.0) nx = nx + 1;
  nz = nx;
  print *,'nx,ny',nx,ny
  do i_Sph = 1,n_Sph 
    x_Sph(i_Sph) = xfrac_Sph(i_Sph)*nx; 
    z_Sph(i_Sph) = zfrac_Sph(i_Sph)*nz; 
    y_Sph(i_Sph) = rad_Sph * (1 - 2. * truncation_by_diameter) + 1
  enddo  

  call system_clock(beginning,rate)
  call Count_and_Allocate
  call make_Lists; print *,'Lists done';
  call Set_om_dependent_variables; 
  call init_u0; 
  it = 0  
  B_accuracy_reached = .false.
  B_dif_ini_found = .false.
  do while (it < nsteps_mx .and. .not.B_accuracy_reached )
    it = it + 1  
    call compute_VelTop; 
    call Stream_and_Collide; 
    if(mod(it,20) == 1) then 
      call calc_dfc;
      if ((it > 500) .and. (abs(dfc_prev - dfc) < abs(dfc_prev_prev - dfc_prev))) then
        if (.not. B_dif_ini_found) then 
          dif_ini = abs(dfc_prev - dfc)  
        endif  
        B_dif_ini_found = .true.
      endif
      print "(I5,F12.3,F12.3,F12.3,F12.3,F12.3,F12.3)",it,real(dfc)       ,aimag(dfc),&
                                                        real(dfc_stress),aimag(dfc_stress),&
                                                        real(dfc_upb),   aimag(dfc_upb)
      if ( (abs(dfc_prev - dfc) / dif_ini < target_ratio) .and. (it > 2000)) then 
        print *,'limit reached',it;  B_accuracy_reached = .true.
      endif  
      dfc_prev_prev = dfc_prev
      dfc_prev = dfc
      open(unit = 20,file='freq_ring_in.asc',action='write',status='unknown', position='append'); 
      write (20,*)   it,real(dfc-dfc_ref),aimag(dfc-dfc_ref); 
      close(unit = 20);   

      open(unit = 20,file='freq_ring_in_stress.asc',action='write',status='unknown', position='append'); 
      write (20,*)   it,real(dfc_stress-dfc_stress_ref),aimag(dfc_stress-dfc_stress_ref); 
      close(unit = 20);   

      open(unit = 20,file='freq_ring_in_upb.asc',action='write',status='unknown', position='append'); 
      write (20,*)   it,real(dfc_upb-dfc_upb_ref),aimag(dfc_upb-dfc_upb_ref); 
      close(unit = 20);   

      open(unit = 20,file='freq_ring_in_noMassLiq.asc',action='write',status='unknown', position='append'); 
      write (20,*)   it,real(dfc_noMassLiq - dfc_noMassLiq_ref),aimag(dfc_noMassLiq - dfc_noMassLiq_ref); 
      close(unit = 20);   

 !     open(unit = 20,file='freq_ring_in_sqr_phase_diss.asc',action='write',status='unknown', position='append'); 
 !     write (20,*)   it,real(dfc), phase_diss
 !     close(unit = 20);   
    endif  
  enddo
  B_save_u = .true.
  call compute_VelTop; call Stream_and_Collide; 
  call write_uxuy
  
  call system_clock(ending); print *,"Elapsed time [secs] : ",real(ending - beginning)/real(rate)
  open(unit = 20,file='sim_time_steps.asc',action='write'); 
  write(20,"(I5,I5,I12,I12)") nx,ny,nPo,nsteps
  write(20,"(F16.8)") real(ending - beginning)/real(rate)
  close(unit = 20);   
  call DeAllocateAll
end program FT_LBM

subroutine init_u0
  use declare; implicit none
  integer :: x,y,z,iPo,iLi
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

logical function in_Sph(x,y,z) 
  use declare; implicit none
  integer :: x,y,z
  in_Sph = .false.
  do i_Sph = 1,n_Sph
    if ( (y.ge.1) .and. ((x-x_Sph(i_Sph))**2 + (y-y_Sph(i_Sph))**2 + (z-z_Sph(i_Sph))**2 .lt. rad_Sph**2) ) then 
      in_Sph = .true.
    endif  
  enddo
end function in_Sph

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

subroutine Set_om_dependent_variables
  use declare; implicit none
  delta_f0 = delta_f0_SI/grid_reso_SI
  f0 = 2.*nu/delta_f0**2 / (2. * pi)
  ovtorder = 2 * iovt - 1; f_SI = ovtorder* f0_SI; 
  delta = delta_f0 / sqrt(float(ovtorder)); 
  omega = 2.*nu/delta**2.; 
  omega2by2 = omega**2/2; Zliq = sqrt(rho0*iu*omega*rho0*nu)
  tauinv_m_omega2by2 = tauinv - omega2by2
  Iom = complex(0.,1.)*omega
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
  allocate( y_arr_nm(ny) )
  allocate( h1(nLi),h2(nLi) )
  allocate( dr(nPo),ux(nPo),uy(nPo),uz(nPo) )
  allocate( uxTop(nPoTop), uzTop(nPoTop) )
  allocate( xBor(nPoBor,0:ndir),yBor(nPoBor,0:ndir),zBor(nPoBor,0:ndir) )
  allocate( Poij(nx,ny,nz),Topj(nx,nz),Linj(nx,ny,nz,0:ndir) )
  allocate( nLiBor(nPoBor),nLi3PLSph(nPo3PL),nLi3PLWal(nPo3PL) )
  allocate( PoLiqx(nPoLiq),PoLiqy(nPoLiq),PoLiqz(nPoLiq) )
  allocate( PoSphx(nPoSph),PoSphy(nPoSph),PoSphz(nPoSph) )
  allocate( PoTopx(nPoTop),PoTopy(nPoTop),PoTopz(nPoTop) )
  allocate( PoBorx(nPoBor),PoBory(nPoBor),PoBorz(nPoBor) )
  allocate( Po3PLx(nPo3PL),Po3PLy(nPo3PL),Po3PLz(nPo3PL) )
  allocate( PoWalLiqx(nPoWalLiq),PoWalLiqy(nPoWalLiq),PoWalLiqz(nPoWalLiq) )
  allocate( permu_Bor(nPoBor,ndir+1),permu_3PL(nPo3PL,ndir+1))
end subroutine Count_and_Allocate

subroutine DeAllocateAll
  use declare; implicit none
  deallocate( y_arr_nm )
  deallocate(h1,h2); 
  deallocate(dr,ux,uy,uz);
  deallocate( uxTop, uzTop )
  deallocate( xBor,yBor,zBor )  
  deallocate( Poij,Linj,Topj )
  deallocate( nLiBor,nLi3PLSph,nLi3PLWal )
  deallocate( PoLiqx,PoLiqy,PoLiqz )
  deallocate( PoSphx,PoSphy,PoSphz )
  deallocate( PoTopx,PoTopy,PoTopz )
  deallocate( PoBorx,PoBory,PoBorz )
  deallocate( Po3PLx,Po3PLy,Po3PLz )
  deallocate( PoWalLiqx,PoWalLiqy,PoWalLiqz )
  deallocate( permu_Bor,permu_3PL)
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
  uxTop_Avg = 0; uzTop_Avg = 0;
  do x = 1,nx
    do z = 1,nz
!      uxTop_Avg = uxTop_Avg - 2*( h1(Linj(x,ny,z,7 )) - h1(Linj(x,ny,z,14)) )/(Zliq-1./3.)
!      uzTop_Avg = uzTop_Avg - 2*( h1(Linj(x,ny,z,11)) - h1(Linj(x,ny,z,17)) )/(Zliq-1./3.)  
       uxTop_Avg = uxTop_Avg + 2*( h1(Linj(x,ny,z,7 )) - h1(Linj(x,ny,z,14)) )/(Zliq+1./3.)
       uzTop_avg = uzTop_Avg + 2*( h1(Linj(x,ny,z,11)) - h1(Linj(x,ny,z,17)) )/(Zliq+1./3.) 
   enddo  
  enddo
  uxTop_Avg = uxTop_Avg/(nx*nz)
  uzTop_Avg = uzTop_Avg/(nx*nz)
  do x = 1,nx
    do z = 1,nz
      uxTop(Topj(x,z)) = uxTop_Avg 
      uzTop(Topj(x,z)) = uzTop_Avg 
    enddo  
  enddo
end subroutine compute_VelTop

subroutine calc_MassDiss
  use declare; implicit none
  integer :: x,y,z,i,iPo
  complex*16 :: ht(0:ndir)
  Mass_Liq = 0;  Qdot_Liq = 0; 
  !$omp parallel private(iPo,i,x,y,z,ht)
    !$omp do
    do iPo = 1,nPoLiq
        x = PoLiqx(iPo)
        y = PoLiqy(iPo)
        z = PoLiqz(iPo)
     do i = 0,ndir
          ht(i) = h1(Linj(x,y,z,i))
     enddo
        call calc_dr_Mass_Diss_loc(ht)
    enddo  
    !$omp end do
  !$omp end parallel

  do iPo = 1,nPoWalLiq
    x = PoWalLiqx(iPo)
    y = PoWalLiqy(iPo)
    z = PoWalLiqz(iPo)
    do i = 0,ndir
      ht(i) = h1(Linj(x,y,z,i))
    enddo
    call calc_dr_Mass_Diss_loc(ht)
  enddo  

  do iPo = 1,nPoBor
    x = PoBorx(iPo)
    y = PoBory(iPo)
    z = PoBorz(iPo)
    do i = 0,ndir
      ht(i) = h1(Linj(x,y,z,i))
    enddo
    call calc_dr_Mass_Diss_loc(ht)
  enddo  

  do iPo = 1,nPo3PL
    x = Po3PLx(iPo)
    y = Po3PLy(iPo)
    z = Po3PLz(iPo)
    do i = 0,ndir
      ht(i) = h1(Linj(x,y,z,i))
    enddo
    call calc_dr_Mass_Diss_loc(ht)
  enddo  

  do iPo = 1,nPoTop
    x = PoTopx(iPo)
    y = PoTopy(iPo)
    z = PoLiqz(iPo)
    do i = 0,ndir
      ht(i) = h1(Linj(x,y,z,i))
    enddo
    call calc_dr_Mass_Diss_loc(ht)
  enddo
end subroutine calc_MassDiss

subroutine calc_dr_Mass_Diss_loc(ht)
  use declare; implicit none
  integer :: i
  complex*16 :: ht(0:ndir),drloc,uxloc,uyloc,uzloc
  complex*16 :: NE_sig_xy, NE_sig_xz, NE_sig_yz, heq, cidotu 
  real*8 :: Qdot_loc 
  uxloc = 0; uyloc = 0; uzloc = 0; 
  do i = 0,ndir
    drloc = drloc + ht(i)
    uxloc = uxloc + ht(i) * dirx(i)
    uyloc = uyloc + ht(i) * diry(i)
    uzloc = uzloc + ht(i) * dirz(i)
  enddo
  Mass_Liq = Mass_Liq + rho0 * abs(uxloc)**2  + abs(uyloc)**2 + abs(uzloc)**2

  NE_sig_xy = 0; NE_sig_xz = 0; NE_sig_yz = 0;
  do i = 0,ndir
    cidotu = dirx(i)*uxloc + diry(i)*uyloc + dirz(i)*uzloc
    heq = wi(i) * (drloc + 3.*cidotu)
    NE_sig_xy = NE_sig_xy - 1./2. * (dirx(i)*diry(i) * (ht(i) - heq))
    NE_sig_xz = NE_sig_xz - 1./2. * (dirx(i)*dirz(i) * (ht(i) - heq))
    NE_sig_yz = NE_sig_yz - 1./2. * (diry(i)*dirz(i) * (ht(i) - heq))
  enddo
  Qdot_loc = 1 / eta0 * (abs(NE_sig_xy)**2 + abs(NE_sig_xz)**2 + abs(NE_sig_yz)**2)
  Qdot_Liq = Qdot_Liq + Qdot_loc
end subroutine calc_dr_Mass_Diss_loc

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

  if (B_save_u) then 
    dr(Poij(x,y,z)) = drloc
    ux(Poij(x,y,z)) = uxloc
    uy(Poij(x,y,z)) = uyloc
    uz(Poij(x,y,z)) = uzloc
  endif  

end subroutine calc_dr_uxyz

subroutine collide(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
  use declare; implicit none
  integer :: x,y,z,i
  complex*16 :: cidotu,heq,drloc,uxloc,uyloc,uzloc,ht(0:ndir)
  cidotu = dirx(i)*uxloc + diry(i)*uyloc + dirz(i)*uzloc
  heq = wi(i) * (drloc + 3.*cidotu)
  h2(Linj(x,y,z,i)) = tauinv_m_omega2by2 * heq - Iom * ht(i)   !omtauinv = 0
end subroutine collide

subroutine Stream_and_Collide
  use omp_lib
  use declare; implicit none
  integer :: x,y,z,i,iloc,ibar,xmd,ymd,zmd,iPo,iLi
  complex*16 :: ht(0:ndir),drloc,uxloc,uyloc,uzloc,cidotu
  FxLiqOnWall = 0; FxLiqOnSphe = 0; 

  !$omp parallel private(iPo,i,x,y,z,drloc,uxloc,uyloc,xmd,ymd,zmd,uzloc,ht)
    !$omp do
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
          call collide(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
        enddo  
    enddo  
    !$omp end do
  !$omp end parallel

  do iPo = 1,nPoWalLiq
    x = PoWalLiqx(iPo)
    y = PoWalLiqy(iPo)
    z = PoWalLiqz(iPo)
    do iloc = 1,5
      ibar = permu_Wal_BBa(iloc)
      i = no_slip(ibar)
      ht(ibar) = h1(Linj(x,y,z,i)) - 6 * wi(i) * dirx(i)  ! page 180 
      FxLiqOnWall = FxLiqOnWall + (ht(ibar) + h1(Linj(x,y,z,i)))*dirx(i)
      
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
      call collide(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
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
      ht(ibar) = h1(Linj(x,y,z,i)) - 6 * wi(i) * cidotu
      FxLiqOnSphe = FxLiqOnSphe + (ht(ibar) + h1(Linj(x,y,z,i)))*dirx(i)
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
      call collide(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
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
      ht(ibar) = h1(Linj(x,y,z,i)) - 6 * wi(i) * cidotu
      FxLiqOnSphe = FxLiqOnSphe + (ht(ibar) + h1(Linj(x,y,z,i)))*dirx(i)
    enddo  
    do iloc = nLi3PLSPh(iPo) + 1, nLi3PLSPh(iPo) + nLi3PLWal(iPo)
      ibar = permu_3PL(iPo,iloc)
      i = no_slip(ibar)
      cidotu = dirx(i)
      ht(ibar) = h1(Linj(x,y,z,i)) - 6 * wi(i) * cidotu
      FxLiqOnWall = FxLiqOnWall + (ht(ibar) + h1(Linj(x,y,z,i)))*dirx(i)
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
      call collide(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
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
      call collide(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
    enddo  
  enddo

  do iLi = 1,nLi
    h1(iLi) = h2(iLi)
  enddo
end subroutine Stream_and_Collide

subroutine write_uxuy
  use declare; implicit none
  integer :: x,y
  open(unit = 20,file='uxuy_3D.asc',action='write')
  write(20,"(I5,I5,I5)")nx,ny,nz
  do x = 1,nx
    do y = 1,ny
      write(20,"(I5,I5,F12.8,F12.8,F12.8,F12.8)") x,y, &
        real(ux(Poij(x,y,int(nz/2.)+1))),aimag(ux(Poij(x,y,int(nz/2.)+1))),&
        real(uy(Poij(x,y,int(nz/2.)+1))),aimag(uy(Poij(x,y,int(nz/2.)+1)))
    enddo
  enddo
  close(unit = 20)
  open(unit = 20,file='ux_vs_y.asc', action='write'); 
  do y = 1,ny
    y_arr_nm(y) = y * grid_reso_SI * 1e9
    write(20,*) y_arr_nm(y),real(ux(Poij(1,y,1))),aimag(ux(Poij(1,y,1)))
  enddo
end subroutine write_uxuy

subroutine calc_dfc
  use declare; implicit none
  real*8 :: df,dGamma, Mass_Wave, Qdot_Wave, Mass_Sph, mq, Zq
  real*8 :: Dt_SI, Dx_SI, u0_SI, Eq_SI, Qdot_SI
  complex*16 :: u_upb,Fx_inertia 
  integer :: x,z
  call calc_MassDiss
  u_upb = 0
  do x = 1,nx
    do z = 1,nz
      u_upb = u_upb + uxTop(Topj(x,z)) 
    enddo  
  enddo
  u_upb = u_upb / nPoTop
 
  Zq = (rho0/rho0_SI)*(nu/nu_SI)*(delta_f0_SI/delta_f0)*Zq_SI
  mq = Zq / (2. * f0); 
  Mass_Sph = n_Sph * rho0 * rad_Sph**3. * 4.*pi/3 * &
      (truncation_by_diameter-1.)**2 * (1 + 2*truncation_by_diameter)  !4/3 (-1 + al)^2 (1 + 2 al) \[Pi]
!  Mass_Sph = nPoSph * rho0

  Mass_Wave = abs(u_upb)**2 * rho0 * delta/2.0  
  Qdot_Wave = abs(u_upb)**2 * aimag(Zliq)

  Dt_SI = f0/f0_SI; Dx_SI = grid_reso_SI; u0_SI = Dx_SI/Dt_SI
  Qdot_SI = (Qdot_Liq/(nx * nz) + Qdot_Wave) * (rho0_SI / rho0) / (nu_SI/nu) * u0_SI**4 * Dx_SI 
  Eq_SI = Zq_SI / (2* f0_SI) * u0_SI**2 / 2. 
  df     = - ovtorder* f0_SI * (Mass_Liq/(nx * nz) + mass_Sph/(nx * nz) + Mass_Wave) / mq
  dGamma = 1./(4.*Pi) * Qdot_SI / Eq_SI
  dfc    = (df + iu * dGamma)
  
  Fx_inertia = -iu*omega*Mass_Sph
  ZL = -(FxLiqOnWall + FxLiqOnSphe + Fx_inertia)/(nx*nz)
  dfc_stress = f0_SI*(iu/pi * ZL/Zq)

  dfc_upb = f0_SI*(iu/pi*Zliq/Zq) * uxTop_Avg / uxTop_Avg_ref
  dfc_noMassLiq = - ovtorder* f0_SI * (mass_Sph/(nx * nz) + Mass_Wave) / mq + iu * dGamma
!  phase_diss = atan2(aimag(Qdot_Liq),real(Qdot_Liq))
end subroutine calc_dfc

