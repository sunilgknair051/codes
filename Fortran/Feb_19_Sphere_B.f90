module declare
  implicit none
  integer, parameter :: ndir = 18
  integer, parameter :: dirx(0:18) = (/0,+1,-1, 0, 0, 0, 0,+1,-1,+1,-1, 0, 0,+1,-1,+1,-1, 0, 0/)
  integer, parameter :: dirz(0:18) = (/0, 0, 0, 0, 0,+1,-1, 0, 0,+1,-1,+1,-1, 0, 0,-1,+1,-1,+1/)
  integer, parameter :: diry(0:18) = (/0, 0, 0,+1,-1, 0, 0,+1,-1, 0, 0,+1,-1,-1,+1, 0, 0,+1,-1/)
  integer, parameter :: no_slip(0:18) = (/0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17/)
  real*8,parameter :: wi(0:18) = (/ 1./3.,1./18,1./18.,1./18.,1./18.,1./18.,1./18., &
                           1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36./)
  real*8,parameter :: pi = 4.D0*DATAN(1.D0), rho0 = 1., nu = 1./6., tauinv_Liq = 2./(6.*nu + 1.) 
  real*8,parameter :: omtauinv_Liq = 1. - tauinv_Liq
  real*8           :: x_Sph, y_Sph, z_Sph, rad_Sph, delta, omega, omega2by2 
  real*8           :: ncycles_mx 
  complex*16 :: Iomega,tauinv_Liq_m_omega2by2
  integer :: nsteps_mx
  integer :: nx,ny,nz,nLi,nPo,nPoLiq,nPoSph,nPoBor
  
  complex*16, allocatable :: dr(:),ux(:),uy(:),uz(:),h1(:), h2(:)

  complex*16, allocatable :: FxBor(:,:),FyBor(:,:),FzBor(:,:)
  real*8    , allocatable :: xBor (:,:), yBor(:,:), zBor(:,:), qBor(:,:)
  complex*16, allocatable :: uxBor(:,:),uyBor(:,:),uzBor(:,:)
  complex*16 :: FxLiqOnSph,ratio
    
  integer, allocatable :: Poij(:,:,:),Linj(:,:,:,:)
  integer, allocatable :: nLiBor(:)
  integer, allocatable :: PoLiqx(:),PoLiqy(:),PoLiqz(:)
  integer, allocatable :: PoSphx(:),PoSphy(:),PoSphz(:)
  integer, allocatable :: PoBorx(:),PoBory(:),PoBorz(:)
  integer, allocatable :: permu_Bor(:,:)
  
  integer :: it,nsteps,error
  logical :: B_Save_uxyz
  integer :: beginning,rate,ending
  
end module declare

program FT_LBM
  use declare; implicit none
  integer :: x,y,z
  
  open(unit = 20,file='freq_ring_in.asc1',action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='uxuy_3D.asc1',     action='write',status='unknown'); close(unit = 20,status='delete');   

  ny = 200; nx = ny; nz = ny; nsteps = 1000;
  delta = int(ny/4.); rad_Sph = ny/10.; x_Sph = nx/2.; z_Sph = nz/2.; y_Sph = ny/2.; 
  omega = 2.*nu / delta**2.; 
  omega2by2 = omega**2 / 2; 
  tauinv_Liq_m_omega2by2 = tauinv_Liq - omega2by2
  Iomega = complex(0.,1.)*omega

  call system_clock(beginning,rate)
  call Count_and_Allocate
  print *,'allocated'
  call make_Lists
  print *,'Lists'
  call init_u0; 
  it = 0  
  B_Save_uxyz = .true.; 
  do while (it < nsteps)
    it = it + 1  
    call Stream_and_Collide; 
    call calc_force_compare
    if (nsteps.gt.10) then 
!      if(mod(it,2) == 1) then 
        print *,it,nsteps,real(ratio),aimag(ratio)
        call write_vels_to_file
!      endif 
!      if (mod(it,5) == 1) then 
        open(unit = 20,file='freq_ring_in.asc1',action='write',status='unknown', position='append'); 
        write (20,*)   it,nsteps,real(ratio),aimag(ratio)
        close(unit = 20);   
!      endif  
    endif  
  enddo
  call system_clock(ending); print *,"Elapsed time is in seconds --->",real(ending - beginning)/real(rate)
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
    ux(Poij(x,y,z)) = 1.
  enddo  
  do iLi = 1,nLi
    h1(iLi) = 0
  enddo
end subroutine init_u0

logical function in_Sph(x,y,z) 
  use declare; implicit none
  integer :: x,y,z
  in_Sph = ((x-x_Sph)**2 + (y-y_Sph)**2 + (z-z_Sph)**2 .lt. rad_Sph**2)
end function in_Sph

logical function in_Liq(x,y,z) 
  use declare; implicit none
  integer :: x,y,z,xmd,ymd,zmd,i
  logical in_Sph,B_mds_not_in_Sph
  B_mds_not_in_Sph = .true.
  do i = 0,ndir
    xmd = mod(nx+x-dirx(i)-1,nx)+1
    ymd = mod(ny+y-diry(i)-1,ny)+1
    zmd = mod(nz+z-dirz(i)-1,nz)+1
    if (in_Sph(xmd,ymd,zmd))  then 
      B_mds_not_in_Sph = .false.
    endif  
  enddo  
  in_Liq = ((.not. in_Sph(x,y,z)) .and. (B_mds_not_in_Sph))
end function in_Liq

logical function on_Bor(x,y,z) 
  use declare; implicit none
  integer :: x,y,z,xmd,ymd,zmd,i
  logical in_Sph,B_Bor
  B_Bor = .false.
  do i = 0,ndir
    xmd = mod(nx+x-dirx(i)-1,nx)+1
    ymd = mod(ny+y-diry(i)-1,ny)+1
    zmd = mod(nz+z-dirz(i)-1,nz)+1
    if ( (.not. in_Sph(x,y,z)) .and. in_Sph(xmd,ymd,zmd) ) then 
      B_Bor = .true.
    endif  
  enddo  
  on_Bor = B_Bor
end function on_Bor

subroutine Count_and_Allocate
  use declare; implicit none
  integer :: x,y,z
  logical :: in_Sph,in_Liq, on_Bor
  nPo = nx*ny*nz; nLi = nPo * (ndir + 1); 
  nPoLiq = 0; nPoSph = 0; nPoBor = 0
  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        if ( in_Liq(x,y,z) ) nPoLiq = nPoLiq + 1        
        if ( in_Sph(x,y,z) ) nPoSph = nPoSph + 1
        if ( on_Bor(x,y,z) ) nPoBor = nPoBor + 1
      enddo
    enddo
  enddo
  allocate( h1(nLi),h2(nLi) )
  allocate( dr(nPo),ux(nPo),uy(nPo),uz(nPo) )
  allocate( FxBor(nPoBor,0:ndir),FyBor(nPoBor,0:ndir),FzBor(nPoBor,0:ndir))
  allocate( qBor(nPoBor,0:ndir), xBor(nPoBor,0:ndir),yBor(nPoBor,0:ndir),zBor(nPoBor,0:ndir) )
  allocate( uxBor(nPoBor,0:ndir),uyBor(nPoBor,0:ndir),uzBor(nPoBor,0:ndir) )
  allocate( Poij(nx,ny,nz),Linj(nx,ny,nz,0:ndir) )
  allocate( nLiBor(nPoBor))
  allocate( PoLiqx(nPoLiq),PoLiqy(nPoLiq),PoLiqz(nPoLiq) )
  allocate( PoSphx(nPoSph),PoSphy(nPoSph),PoSphz(nPoSph) )
  allocate( PoBorx(nPoBor),PoBory(nPoBor),PoBorz(nPoBor) )
  allocate( permu_Bor(ndir+1,nPoBor))
end subroutine Count_and_Allocate

subroutine DeAllocateAll
  use declare; implicit none
  deallocate( h1,h2 ); 
  deallocate( dr,ux,uy,uz );
  deallocate( FxBor,FyBor,FzBor )
  deallocate( qBor,xBor,yBor,zBor )  
  deallocate( uxBor,uyBor,uzBor )  
  deallocate( Poij,Linj )
  deallocate( nLiBor )
  deallocate( PoLiqx,PoLiqy,PoLiqz )
  deallocate( PoSphx,PoSphy,PoSphz )
  deallocate( PoBorx,PoBory,PoBorz )
  deallocate( permu_Bor)
end subroutine DeAllocateAll

subroutine make_Lists
  use declare; implicit none
  integer :: x,y,z,xmd,ymd,zmd,i,i_next
  logical :: in_Sph, in_Liq, on_Bor
  nPo = 0; nLi = 0; nPoLiq = 0; nPoSph = 0; 
  nPoBor = 0; 
  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        nPo = nPo + 1; 
        Poij(x,y,z) = nPo
        do i = 0,ndir
          nLi = nLi + 1 
          Linj(x,y,z,i) = nLi
        enddo

        if ( in_Liq(x,y,z) )  then
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

        if (on_Bor(x,y,z)) then
          nPoBor = nPoBor + 1
          PoBorx(nPoBor) = x
          PoBory(nPoBor) = y
          PoBorz(nPoBor) = z
          nLiBor(nPoBor) = 0  
          do i = 0,ndir
            xmd = mod(nx+x-dirx(i)-1,nx)+1
            ymd = mod(ny+y-diry(i)-1,ny)+1
            zmd = mod(nz+z-dirz(i)-1,nz)+1
            if (in_Sph(xmd,ymd,zmd)) then 
              nLiBor(nPoBor) = nLiBor(nPoBor) + 1
              permu_Bor(nLiBor(nPoBor),nPoBor) = i
            endif
          enddo
          i_next = nLiBor(nPoBor)
          do i = 0,ndir
            xmd = mod(nx+x-dirx(i)-1,nx)+1
            ymd = mod(ny+y-diry(i)-1,ny)+1
            zmd = mod(nz+z-dirz(i)-1,nz)+1
            if ( .not.in_Sph(xmd,ymd,zmd) )  then 
              i_next = i_next + 1
              permu_Bor(i_next,nPoBor) = i
            endif
          enddo
        endif  
      enddo
    enddo
  enddo
end subroutine make_Lists

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
  if (B_Save_uxyz) then
    dr(Poij(x,y,z)) = drloc
    ux(Poij(x,y,z)) = uxloc
    uy(Poij(x,y,z)) = uyloc
    uz(Poij(x,y,z)) = uzloc
  endif
end subroutine calc_dr_uxyz

subroutine Stream_and_Collide
  use omp_lib
  use declare; implicit none
  integer :: x,y,z,i,iloc,ibar,xmd,ymd,zmd,iPo,iLi
  complex*16 :: ht(0:ndir),drloc,uxloc,uyloc,uzloc,cidotu
  do iPo = 1,nPoBor
    do iloc = 0,ndir
      FxBor(iPo,iloc) = 0
      FyBor(iPo,iloc) = 0
      FzBor(iPo,iloc) = 0
    enddo
  enddo  

!  !$omp parallel private(iPo,i,x,y,z,drloc,uxloc,uyloc,xmd,ymd,zmd,uzloc,ht)
!    !$omp do
      do iPo = 1,nPoLiq
        x = PoLiqx(iPo)
        y = PoLiqy(iPo)
        z = PoLiqz(iPo)
        do i = 0,ndir
          xmd = mod(nx+x-dirx(i)-1,nx)+1
          ymd = mod(ny+y-diry(i)-1,ny)+1
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
  FxLiqOnSph = 0
  do iPo = 1,nPoBor
    x = PoBorx(iPo)
    y = PoBory(iPo)
    z = PoBorz(iPo)
    do iloc = 1,nLiBor(iPo)
      i = permu_Bor(iloc,iPo)
      xmd = mod(nx+x-dirx(i)-1,nx)+1
      ymd = mod(ny+y-diry(i)-1,ny)+1
      zmd = mod(nz+z-dirz(i)-1,nz)+1
      ibar = no_slip(i)
      cidotu = dirx(ibar)
      ht(i) = h1(Linj(x,y,z,ibar)) - 6 * wi(ibar) * cidotu
      FxBor(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * dirx(i)
      FxLiqOnSph = FxLiqOnSph + FxBor(iPo,iloc)
    enddo  
    do iloc = nLiBor(iPo)+1,ndir+1
      i = permu_Bor(iloc,iPo)
      xmd = mod(nx+x-dirx(i)-1,nx)+1
      ymd = mod(ny+y-diry(i)-1,ny)+1
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

subroutine calc_force_compare
  use declare; implicit none
  real*8 :: FOS_analy_Re,FOS_analy_Im, eta
  eta = nu * rho0
  FOS_analy_Re = 6. * pi * eta * rad_Sph *(1. + rad_Sph/delta)
  FOS_analy_Im = 4.*pi/3.*rho0 * rad_Sph**3 * (1./2. + 9./4.*delta/rad_Sph)
  ratio = (FOS_analy_Re + complex(0,1)*omega * FOS_analy_Im) / FxLiqOnSph
end subroutine calc_force_compare

subroutine write_vels_to_file
  use declare; implicit none
  integer :: x,y,z
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
end subroutine write_vels_to_file
