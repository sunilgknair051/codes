module declare
  implicit none
  real*8, parameter :: pi = 4.D0*DATAN(1.D0)
  real*8, parameter :: f0_SI = 5e6, nu_SI = 1e-6, delta_f0_SI = (2*nu_SI/(2*pi*f0_SI))**0.5
  real*8,parameter :: Zq_SI = 8.8e6, rho0_SI = 1e3, df_analy = 716.799 
  real*8, parameter :: nu = 1./6., tau = 3. * nu + 0.5, rho0 = 1
  real*8, parameter :: tauinv = 2.0/(6.0*nu+1.0), omtauinv = 1.0-tauinv
  real*8 :: cell_height_SI, grid_reso_SI, rad_Sph_SI, rho_Sph_SI
  real*8 :: omega, delta_f0, delta
  integer, parameter :: n_Sph_max = 10
  integer :: i_Sph,n_Sph
  real*8 :: x_Sph(n_Sph_max), y_Sph(n_Sph_max), z_Sph(n_Sph_max), xfrac_Sph(n_Sph_max), zfrac_Sph(n_Sph_max) 
  real*8 :: rad_Sph, nufac_Re_Sph, nufac_Im_Sph
  real*8 :: truncation_by_diameter, cell_diam_by_sphere_diam

  real*8 :: ncycles_mx, target_acc_Hz
  real*8 :: Lambda 
  real*8 :: om_sym, om_sym_m_omega2by2
  real*8 :: om_asy, om_asy_m_omega2by2
  complex*16 :: om_sym_Sph, om_sym_m_omega2by2_Sph, unity_m_om_sym_m_Iom_Sph, unity_m_om_sym_m_Iom  
  complex*16 :: om_asy_Sph, om_asy_m_omega2by2_Sph, unity_m_om_asy_m_Iom_Sph, unity_m_om_asy_m_Iom  
  
  complex*16 :: nuc_Sph, tauinv_Sph, omtauinv_Sph, tauinv_Sph_m_omega2by2
  complex*16 :: tauinv_m_omega2by2,omtauinv_m_Iom,omtauinv_Sph_m_Iom 
  complex*16 :: dfc, dfc_ref, dfc_dif, Iom, Zliq 
  integer :: nx, ny, nz, iovt, ovtorder
  integer :: it, nsteps_mx
  integer, parameter :: ndir = 18
  integer, parameter :: dirx(0:18) = (/0,+1,-1, 0, 0, 0, 0,+1,-1,+1,-1, 0, 0,+1,-1,+1,-1, 0, 0/)
  integer, parameter :: dirz(0:18) = (/0, 0, 0, 0, 0,+1,-1, 0, 0,+1,-1,+1,-1, 0, 0,-1,+1,-1,+1/)
  integer, parameter :: diry(0:18) = (/0, 0, 0,+1,-1, 0, 0,+1,-1, 0, 0,+1,-1,-1,+1, 0, 0,+1,-1/)
  integer, parameter :: no_slip(0:18) = (/0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17/)
  real*8,  parameter :: wi(0:18) = (/ 1./3.,1./18,1./18.,1./18.,1./18.,1./18.,1./18., &
                           1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36./)
  integer :: nLi,nPo,nPoLiq,nPoWalLiq,nPoTop,nPoWalSph,nPoSph
  integer :: beginning,rate,ending
  logical :: B_Save_uxyz

  complex*16, allocatable :: dr(:),ux(:),uy(:),uz(:),h1(:),h2(:)
  complex*16, allocatable :: uxTop(:), uzTop(:)
  complex*16, allocatable :: FxOnWalLiq(:,:),FxOnWalSph(:,:)
  real*8, allocatable :: y_arr_nm(:)
    
  integer, allocatable :: Poij(:,:,:),Linj(:,:,:,:),Topj(:,:)
  integer, allocatable :: PoLiqx(:),PoLiqy(:),PoLiqz(:)
  integer, allocatable :: PoSphx(:),PoSphy(:),PoSphz(:)
  integer, allocatable :: PoTopx(:),PoTopy(:),PoTopz(:)
  integer, allocatable :: PoWalLiqx(:),PoWalLiqy(:),PoWalLiqz(:)
  integer, allocatable :: PoWalSphx(:),PoWalSphy(:),PoWalSphz(:)
  integer, parameter :: permu_Top_BBa(1:5 ) = (/        4,      8,        13,12,            18/) !page87
  integer, parameter :: permu_Top_Liq(1:14) = (/0,1,2,3,  5,6,7,  9,10,11,      14,15,16,17   /)
  integer, parameter :: permu_Wal_BBA(1:5 ) = (/      3,      7,       14,      11,      17   /)
  integer, parameter :: permu_Wal_Liq(1:14) = (/0,1,2,  4,5,6,  8,9,10,   12,13,   15,16,   18/)
end module declare

program FT_LBM
  use declare; implicit none
  complex * 16 :: dfc_prev 
  logical :: B_accuracy_reached
  
  open(unit = 20,file='freq_ring_in.asc',action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='ux_vs_y.asc',     action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='uxuy_3D.asc',     action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='sim_pars.asc',    action='write',status='unknown'); close(unit = 20,status='delete');   
  
  open(unit = 21,file='sim_config.asc', action='read',status='unknown'); 
  read(21,*) target_acc_Hz
  read(21,*) rad_Sph_SI
  read(21,*) cell_height_SI
  read(21,*) grid_reso_SI
  read(21,*) cell_diam_by_sphere_diam
  read(21,*) truncation_by_diameter
  read(21,*) nufac_Re_Sph
  read(21,*) nufac_Im_Sph
  read(21,*) Lambda
  read(21,*) iovt
  read(21,*) n_Sph
  do i_Sph = 1,n_Sph 
    read(21,*) xfrac_Sph(i_Sph),zfrac_Sph(i_Sph)
  enddo
  close(unit = 21);   

  delta_f0 = delta_f0_SI / grid_reso_SI; 
  B_Save_uxyz = .false.;!  target_acc_Hz = 1;

!  rad_Sph_SI = 400e-9;
!  cell_height_SI = 1200e-9; 
!  grid_reso_SI = 40e-9;
!  cell_diam_by_sphere_diam = 1.5;
!  truncation_by_diameter = 0.5;
!  nufac_Re_Sph = 1.;
!  nufac_Im_Sph = 2.;
!  Lambda = 0.25;

  nuc_Sph = complex(nufac_Re_Sph * nu, nufac_Im_Sph * nu)
  tauinv_Sph = 2.0/(6.0*nuc_Sph + 1.0)
  omtauinv_Sph = 1.0 - tauinv_Sph
  ncycles_mx = 2
  iovt = 1
  ny = int(cell_height_SI/grid_reso_SI)
  rad_Sph = 0; nx = 1; nz = nx; 
  call Count_and_Allocate
  call make_Lists
  call Set_om_dependent_variables
!  call init_u0; 
!  do it = 1,nsteps_mx
!    call Stream_and_Collide; 
!    call calc_dfc;
!    if(mod(it,int(nsteps_mx/5)) == 1) then 
!      print "(I9,F15.5,F15.5)", it,real(dfc),aimag(dfc)
!    endif 
!  enddo
  dfc_ref = dfc;
  call DeAllocateAll

  rad_Sph = rad_Sph_SI / grid_reso_SI
  nx = int(2. * rad_Sph * cell_diam_by_sphere_diam); 
  if (mod(nx,2) == 0) nx = nx + 1; nz = nx; 
  do i_Sph = 1,n_Sph 
    x_Sph(i_Sph) = xfrac_Sph(i_Sph)*nx; 
    z_Sph(i_Sph) = zfrac_Sph(i_Sph)*nz; 
    y_Sph(i_Sph) = rad_Sph * (1 - 2. * truncation_by_diameter) + 1
  enddo  

  call system_clock(beginning,rate)
  call Count_and_Allocate
  call make_Lists
  call Set_om_dependent_variables; 
  call init_u0; 
  it = 0  
  B_accuracy_reached = .false.
  do while (it < nsteps_mx .and. .not.B_accuracy_reached )
    it = it + 1  
    call Stream_and_Collide; 
    if(mod(it,100) == 1) then 
      call calc_dfc;
      print "(I5,I7,F16.5,F16.5)",it,nsteps_mx,real(dfc),aimag(dfc)
      if ( abs(dfc_prev - dfc) < target_acc_Hz ) then 
        print *,'limit reached',it;  B_accuracy_reached = .true.
      endif  
      dfc_prev = dfc
    endif  
    if(mod(it,10) == 1) then 
      call calc_dfc;
      open(unit = 20,file='freq_ring_in.asc',action='write',status='unknown', position='append'); 
      write (20,*)   it,real(dfc-dfc_ref),aimag(dfc-dfc_ref); 
      close(unit = 20);   
    endif  
  enddo
  B_Save_uxyz = .true.; call Stream_and_Collide; B_Save_uxyz = .false.; 
  call system_clock(ending); print *,"Elapsed time [secs]: ",real(ending - beginning)/real(rate)
  open(unit = 20,file='sim_time_steps.asc',action='write'); 
  write(20,"(I5,I5,I12,I12)") nx,ny,nPo,it
  write(20,"(F12.8)") real(ending - beginning)/real(rate)
  close(unit = 20);   
  call Plot_vels;
  call DeAllocateAll
end program FT_LBM

!subroutine Set_pars_for_TRT
!  use declare; implicit none
!
!  om_sym = 1 / (3. * nu + 1./2.) 
!  om_asy = (Lambda / (1./om_sym - 1./2.) + 1./2.)**(-1)
!  om_sym_m_omega2by2 = om_sym - omega**2/2.
!  om_asy_m_omega2by2 = om_asy - omega**2/2.
!  unity_m_om_sym_m_Iom = 1 - om_sym - Iom 
!  unity_m_om_asy_m_Iom = 1 - om_asy - Iom 
!
!  om_sym_Sph = 1 / (3. * nuc_Sph + 1./2.) 
!  om_asy_Sph = (Lambda / (1./om_sym_Sph - 1./2.) + 1./2.)**(-1)
!  om_sym_m_omega2by2_Sph = om_sym_Sph - omega**2/2.
!  om_asy_m_omega2by2_Sph = om_asy_Sph - omega**2/2.
!  unity_m_om_sym_m_Iom_Sph = 1 - om_sym_Sph - Iom 
!  unity_m_om_asy_m_Iom_Sph = 1 - om_asy_Sph - Iom 
!end subroutine Set_pars_for_TRT

subroutine Set_om_dependent_variables
  use declare; implicit none
  real*8 :: eta_Sph
  ovtorder = 2 * iovt - 1; 
  delta = delta_f0 / sqrt(float(ovtorder)); 
  omega = 2.*nu/delta**2.; 
  Iom = complex(0,1)*omega
  omtauinv_m_Iom     = omtauinv     - Iom
  omtauinv_Sph_m_Iom = omtauinv_Sph - Iom
  Zliq = -sqrt(rho0*Iom*rho0*nu)
  tauinv_m_omega2by2     = tauinv     - omega**2/2.
  tauinv_Sph_m_omega2by2 = tauinv_Sph - omega**2/2.
  eta_Sph = nufac_Re_Sph * nu * rho0
  nsteps_mx = int(ncycles_mx / omega * (ny / delta)**2)
end subroutine Set_om_dependent_variables

subroutine init_u0
  use declare; implicit none
  integer :: iLi
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
      in_Sph = .true.
    endif  
  enddo
end function in_Sph

subroutine Count_and_Allocate
  use declare; implicit none
  integer :: x,y,z
  logical :: in_Sph
  nPo = nx*ny*nz; nLi = nPo * (ndir + 1); 
  nPoTop = nx*nz; nPoWalLiq = 0; 
  nPoWalSph = 0;
  nPoLiq = 0;  nPoSph = 0; 
  do x = 1,nx
    do y = 2,ny-1 !!!
      do z = 1,nz
        if (in_Sph(x,y,z)) then 
          nPoSph = nPoSph + 1                
        else 
          nPoLiq = nPoLiq + 1        
        endif
      enddo  
    enddo
  enddo
  do x = 1,nx
    do z = 1,nz
      if (in_Sph(x,y,z)) then 
        nPoWalSph = nPoWalSph + 1
      else 
        nPoWalLiq = nPoWalLiq + 1        
      endif  
    enddo  
  enddo

  allocate( y_arr_nm(ny) )
  allocate( h1(nLi),h2(nLi) )
  allocate( dr(nPo),ux(nPo),uy(nPo),uz(nPo) )
  allocate( uxTop(nPoTop), uzTop(nPoTop) )
  allocate( Poij(nx,ny,nz),Topj(nx,nz),Linj(nx,ny,nz,0:ndir) )
  allocate( FxOnWalLiq(nPoWalLiq,0:ndir), FxOnWalSph(nPoWalSph,0:ndir)  )
  allocate( PoLiqx(nPoLiq),PoLiqy(nPoLiq),PoLiqz(nPoLiq) )
  allocate( PoSphx(nPoSph),PoSphy(nPoSph),PoSphz(nPoSph) )
  allocate( PoTopx(nPoTop),PoTopy(nPoTop),PoTopz(nPoTop) )
  allocate( PoWalLiqx(nPoWalLiq),PoWalLiqy(nPoWalLiq),PoWalLiqz(nPoWalLiq) )
  allocate( PoWalSphx(nPoWalSph),PoWalSphy(nPoWalSph),PoWalSphz(nPoWalSph) )
end subroutine Count_and_Allocate

subroutine DeAllocateAll
  use declare; implicit none
  deallocate( y_arr_nm )
  deallocate( h1,h2 ); 
  deallocate( dr,ux,uy,uz );
  deallocate( uxTop, uzTop )
  deallocate( Poij,Topj,Linj )
  deallocate( FxOnWalLiq, FxOnWalSph )
  deallocate( PoLiqx,PoLiqy,PoLiqz )
  deallocate( PoSphx,PoSphy,PoSphz )
  deallocate( PoTopx,PoTopy,PoTopz )
  deallocate( PoWalLiqx,PoWalLiqy,PoWalLiqz )
  deallocate( PoWalSphx,PoWalSphy,PoWalSphz )
end subroutine DeAllocateAll


subroutine make_Lists
  use declare; implicit none
  integer :: x,y,z,i
  logical in_Sph
  nPo = 0; nLi = 0; nPoLiq = 0; nPoSph = 0; nPoTop = 0; nPoWalLiq = 0; nPoWalSph = 0; 
  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        nPo = nPo + 1; 
        Poij(x,y,z) = nPo
        do i = 0,ndir
          nLi = nLi + 1 
          Linj(x,y,z,i) = nLi
        enddo
      enddo
    enddo
  enddo  
  do x = 1,nx
    do y = 2,ny-1  !!!
      do z = 1,nz
        if (.not.in_Sph(x,y,z))  then 
          nPoLiq = nPoLiq + 1        
          PoLiqx(nPoLiq) = x
          PoLiqy(nPoLiq) = y
          PoLiqz(nPoLiq) = z
        endif   
        if (in_Sph(x,y,z)) then 
          nPoSph = nPoSph + 1
          PoSphx(nPoSph) = x
          PoSphy(nPoSph) = y
          PoSphz(nPoSph) = z
        endif  
      enddo
    enddo
  enddo  

  do x = 1,nx
    do z = 1,nz
      if (.not.in_Sph(x,y,z)) then 
        nPoWalLiq = nPoWalLiq + 1
        PoWalLiqx(nPoWalLiq) = x
        PoWalLiqy(nPoWalLiq) = 1
        PoWalLiqz(nPoWalLiq) = z
      endif  
        if (in_Sph(x,y,z)) then 
        nPoWalSph = nPoWalSph + 1
        PoWalSphx(nPoWalSph) = x
        PoWalSphy(nPoWalSph) = 1
        PoWalSphz(nPoWalSph) = z
      endif  
      nPoTop = nPoTop + 1
      Topj(x,z) = nPoTop
      PoTopx(nPoTop) = x
      PoTopy(nPoTop) = ny
      PoTopz(nPoTop) = z
    enddo
  enddo  
end subroutine make_Lists

subroutine calc_dr_uxyz(x,y,z,ht,drloc,uxloc,uyloc,uzloc)
  use declare; implicit none
  complex*16 :: drloc,ht(0:ndir),uxloc,uyloc,uzloc
  integer :: x,y,z,i
  drloc = 0; uxloc = 0; uyloc = 0; uzloc = 0; 
  do i = 0,ndir
    drloc = drloc + ht(i)
    uxloc = uxloc + ht(i) * dirx(i)
    uyloc = uyloc + ht(i) * diry(i)
    uzloc = uzloc + ht(i) * dirz(i)
  enddo  
  if (B_Save_uxyz) then
    dr(Poij(x,y,z)) = drloc
    ux(Poij(x,y,z)) = uxloc
    uy(Poij(x,y,z)) = uyloc
    uz(Poij(x,y,z)) = uzloc
  endif
end subroutine calc_dr_uxyz

subroutine collide(x,y,z,i,ht,drloc,uxloc,uyloc,uzloc,tauinv_m_omega2by2_loc,omtauinv_m_Iom_loc)
  use declare; implicit none
  integer :: x,y,z,i
  complex*16 :: cidotu,heq,drloc,uxloc,uyloc,uzloc,ht(0:ndir)
  complex*16 :: tauinv_m_omega2by2_loc,omtauinv_m_Iom_loc
  cidotu = dirx(i)*uxloc + diry(i)*uyloc + dirz(i)*uzloc
  heq = wi(i) * (drloc + 3.0*cidotu)
  h2(Linj(x,y,z,i)) = tauinv_m_omega2by2_loc * heq + omtauinv_m_Iom_loc * ht(i)   !omtauinv = 0
end subroutine collide

subroutine collide_TRT(x,y,z,ht,drloc,uxloc,uyloc,uzloc,&
    om_sym_m_omega2by2_loc, unity_m_om_sym_m_Iom_loc,&
    om_asy_m_omega2by2_loc, unity_m_om_asy_m_Iom_loc)
  use declare; implicit none
  integer,parameter :: n_sa = 9
  integer :: x,y,z,i,j
  complex*16 :: cidotu,drloc,uxloc,uyloc,uzloc,ht(0:ndir),heqi(0:ndir)
  complex*16 :: om_sym_m_omega2by2_loc,unity_m_om_sym_m_Iom_loc
  complex*16 :: om_asy_m_omega2by2_loc,unity_m_om_asy_m_Iom_loc
  complex*16 :: h_sym(0:9),h_asy(0:9),heq_sym(0:9),heq_asy(0:9),ht_sym(0:9),ht_asy(0:9)
  integer, parameter :: j2i(0:9   ) = (/0,1,  3,  5,  7,  9,  11, 13, 15, 17/)
  do i = 0,ndir                  
    cidotu = dirx(i)*uxloc + diry(i)*uyloc + dirz(i)*uzloc
    heqi(i) = wi(i) * (drloc + 3.0*cidotu)
  enddo  
  do j = 0,n_sa                  
    h_sym(    j) = ( ht(  j2i(j)) + ht(  no_slip(j2i(j))) / 2. )
    h_sym(    j) = ( ht(  j2i(j)) - ht(  no_slip(j2i(j))) / 2. )
    heq_sym(  j) = ( heqi(j2i(j)) + heqi(no_slip(j2i(j))) / 2. )
    heq_asy(  j) = ( heqi(j2i(j)) - heqi(no_slip(j2i(j))) / 2. )
    ht_sym(j) = om_sym_m_omega2by2_loc * heq_sym(j) + unity_m_om_sym_m_Iom_loc * h_sym(j)
    ht_asy(j) = om_asy_m_omega2by2_loc * heq_asy(j) + unity_m_om_asy_m_Iom_loc * h_asy(j)
    ht_sym(j) = unity_m_om_sym_m_Iom * h_sym(j) + om_sym_m_omega2by2 * heq_sym(j)
    ht_asy(j) = unity_m_om_asy_m_Iom * h_asy(j) + om_asy_m_omega2by2 * heq_asy(j) 
  enddo  
  do j = 0,n_sa                  
    h2(Linj(x,y,z,        j2i(j) )) = ht_sym(j) + ht_asy(j)  
    h2(Linj(x,y,z,no_slip(j2i(j)))) = ht_sym(j) - ht_asy(j)   
  enddo  
end subroutine collide_TRT


!subroutine collide_TRT(x,y,z,ht,drloc,uxloc,uyloc,uzloc,&
!    om_sym_m_omega2by2_loc, unity_m_om_sym_m_Iom_loc,&
!    om_asy_m_omega2by2_loc, unity_m_om_asy_m_Iom_loc)
!  use declare; implicit none
!  integer :: x,y,z,i,j
!  complex*16 :: cidotu,drloc,uxloc,uyloc,uzloc,ht(0:ndir),heqi(0:ndir)
!  complex*16 :: om_sym_m_omega2by2_loc,unity_m_om_sym_m_Iom_loc
!  complex*16 :: om_asy_m_omega2by2_loc,unity_m_om_asy_m_Iom_loc
!  integer, parameter :: n_sa = 9 
!  integer, parameter :: j2i(0:9   ) = (/0,1,  3,  5,  7,  9,  11, 13, 15, 17/)
!  integer, parameter :: i2j(0:ndir) = (/0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9/)
!  complex*16 :: h_sym(0:9),h_asy(0:9),heq_sym(0:9),heq_asy(0:9),ht_sym(0:9),ht_asy(0:9)
!  do i = 0,ndir                  
!    cidotu = dirx(i)*uxloc + diry(i)*uyloc + dirz(i)*uzloc
!    heqi(i) = wi(i) * (drloc + 3.0*cidotu)
!  enddo  
!  do j = 0,n_sa                  
!    h_sym(  j) = ( ht(  j) + ht(  j2i(j)) / 2. )
!    h_asy(  j) = ( ht(  j) - ht(  j2i(j)) / 2. )
!    heq_sym(j) = ( heqi(j) + heqi(j2i(j)) / 2. )
!    heq_asy(j) = ( heqi(j) - heqi(j2i(j)) / 2. )
!    ht_sym(j) = om_sym_m_omega2by2_loc * heq_sym(j) + unity_m_om_sym_m_Iom_loc * h_sym(j)
!    ht_asy(j) = om_asy_m_omega2by2_loc * heq_asy(j) + unity_m_om_asy_m_Iom_loc * h_asy(j)
!  enddo  
!  h2(Linj(x,y,z,0)) = ht_sym(0)
!  do j = 1,n_sa                  
!    h2(Linj(x,y,z,j  )) = ( ht_sym(j2i(j)) + ht_asy(j2i(j)) ) 
!    h2(Linj(x,y,z,j+1)) = ( ht_sym(j2i(j)) - ht_asy(j2i(j)) )   
!  enddo  
!end subroutine collide_TRT
!

subroutine Stream_and_Collide
  use omp_lib
  use declare; implicit none
  integer :: x,y,z,i,iloc,ibar,xmd,ymd,zmd,iPo,iLi
  complex*16 :: ht(0:ndir),drloc,uxloc,uyloc,uzloc,cidotu

  !!$omp parallel private(iPo,x,y,z,i,drloc,uxloc,uyloc,uzloc,xmd,ymd,zmd,ht)
    !!$omp do
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
          call collide(x,y,z,i,ht,drloc,uxloc,uyloc,uzloc,tauinv_m_omega2by2,omtauinv_m_Iom)
!          call collide_TRT(x,y,z,ht,drloc,uxloc,uyloc,uzloc,&
!                om_sym_m_omega2by2, unity_m_om_sym_m_Iom,&
!                om_asy_m_omega2by2, unity_m_om_asy_m_Iom)          
        enddo  
      enddo  
    !!$omp end do  
  !!$omp end parallel
  !!$omp parallel private(iPo,x,y,z,i,drloc,uxloc,uyloc,uzloc,xmd,ymd,zmd,ht)
    !!$omp do
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
          call collide(x,y,z,i,ht,drloc,uxloc,uyloc,uzloc,tauinv_Sph_m_omega2by2,omtauinv_Sph_m_Iom)
        enddo  
      enddo  
    !!$omp end do  
  !!$omp end parallel

  do iPo = 1,nPoWalLiq
    x = PoWalLiqx(iPo)
    y = PoWalLiqy(iPo)
    z = PoWalLiqz(iPo)
    do iloc = 1,5
      ibar = permu_Wal_BBa(iloc)
      i = no_slip(ibar)
      ht(ibar) = h1(Linj(x,y,z,i)) - 6 * wi(i) * dirx(i)  ! page 180 
      FxOnWalLiq(iPo,iloc) = (ht(ibar) + h1(Linj(x,y,z,i))) * dirx(i)   ! page 217
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
      call collide(x,y,z,i,ht,drloc,uxloc,uyloc,uzloc,tauinv_m_omega2by2,omtauinv_m_Iom)
    enddo  
  enddo  

  do iPo = 1,nPoWalSph
    x = PoWalSphx(iPo)
    y = PoWalSphy(iPo)
    z = PoWalSphz(iPo)
    do iloc = 1,5
      ibar = permu_Wal_BBa(iloc)
      i = no_slip(ibar)
      ht(ibar) = h1(Linj(x,y,z,i)) - 6 * wi(i) * dirx(i)  ! page 180 
      FxOnWalSph(iPo,iloc) = (ht(ibar) + h1(Linj(x,y,z,i))) * dirx(i)   ! page 217
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
      call collide(x,y,z,i,ht,drloc,uxloc,uyloc,uzloc,tauinv_Sph_m_omega2by2,omtauinv_Sph_m_Iom)
    enddo  
  enddo  

  do iPo = 1,nPoTop
    x = PoTopx(iPo)
    y = PoTopy(iPo)
    z = PoTopz(iPo)
    do iloc = 1,5
      i = permu_Top_BBa(iloc)
      ibar = no_slip(i)
      uxTop(Topj(x,z)) = -2*wi(7)*( h1(Linj(x,ny,z,7 )) - h1(Linj(x,ny,z,14)) )/(Zliq-1./3.)
      cidotu = dirx(ibar)*uxTop(Topj(x,z)) 
      ht(i) = h1(Linj(x,y,z,ibar)) - 6 * wi(i) * cidotu
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
      call collide(x,y,z,i,ht,drloc,uxloc,uyloc,uzloc,tauinv_m_omega2by2,omtauinv_m_Iom)
    enddo  
  enddo
  do iLi = 1,nLi
    h1(iLi) = h2(iLi)
  enddo
end subroutine Stream_and_Collide

subroutine calc_dfc
  use declare; implicit none
  integer :: iPo,iloc
  complex*16 :: FxTotOnWal,ZL
  real *8 :: Zq
  FxTotOnWal = 0.
  do iPo = 1,nPoWalLiq
    do iloc = 1,5
      FxTotOnWal = FxTotOnWal + FxOnWalLiq(iPo,iloc)     
    enddo
  enddo
  do iPo = 1,nPoWalSph
    do iloc = 1,5
      FxTotOnWal = FxTotOnWal + FxOnWalSph(iPo,iloc)     
    enddo
  enddo
  ZL = (FxTotOnWal)/(nx*nz)
  Zq = (rho0/rho0_SI)*(nu/nu_SI)*(delta_f0_SI/delta_f0)*Zq_SI
  dfc    = - f0_SI * complex(0,1)/(pi*Zq) * ZL
end subroutine calc_dfc

subroutine Plot_vels
  use declare; implicit none
  integer :: x,y
  open(unit = 20,file='uxuy_3D.asc',action='write')
  write(20,"(I5,I5)")nx,ny
  do x = 1,nx
    do y = 1,ny
      write(20,"(I5,I5,F17.8,F17.8,F17.8,F17.8)") x,y, &
                                 real(ux(Poij(x,y,int(nz/2)))),aimag(ux(Poij(x,y,int(nz/2)))),&
                                 real(uy(Poij(x,y,int(nz/2)))),aimag(uy(Poij(x,y,int(nz/2)))) 
    enddo
  enddo
  close(unit = 20)
  open(unit = 20,file='ux_vs_y.asc', action='write'); 
  do y = 1,ny
    y_arr_nm(y) = y * grid_reso_SI * 1e9
    write(20,*) y_arr_nm(y),real(ux(Poij(1,y,1))),aimag(ux(Poij(1,y,1)))
  enddo
end subroutine Plot_vels
