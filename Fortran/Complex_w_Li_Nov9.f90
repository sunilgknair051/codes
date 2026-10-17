module declare
  implicit none
  integer,parameter :: nx=2, ny=50, ndir=8, nPo=nx*ny, nPoLi=nx*ny*(ndir+1)
  integer,parameter :: dirx(0:8)=(/0,1,0,-1,0,1,-1,-1,1/)
  integer,parameter :: diry(0:8)=(/0,0,1,0,-1,1,1,-1,-1/)
  integer,parameter :: no_slip(0:8)=(/0,3,4,1,2,7,8,5,6/)
  real,parameter :: w0 = 4./9., ws=1./9., wd= 1./36.
  real,parameter :: wi(0:8)=(/w0,ws,ws,ws,ws,wd,wd,wd,wd/)
  real,parameter :: nu = 1./6., rho0 = 1., tauinv = 2./(6.*nu + 1.), omtauinv = 1. - tauinv
  real,parameter :: x_sphere =nx/2., y_sphere = -1., rad_sphere = 0.
  real,parameter :: PI = 4.D0*DATAN(1.D0)
  real,parameter :: delfac = 1./6., delta = ny*delfac, omega = 2.*nu/delta**2., omega2by2 = omega**2./2.
  integer,parameter :: NSTEPS = 105
  ! integer,parameter :: NSTEPS = int(2 / omega)
  real,parameter :: uxWall = 1., Zliq = rho0 *(omega*nu)**(1/2)
  real,parameter :: Zq_SI = 8.8E6, nu_SI = 1E-6, rho_Liq_SI = 1e3, f0_SI = 5E6, delta_SI = 250E-9
  real,parameter :: dx_SI = 250e-9/delta, omega_SI = 2. * Pi *f0_SI
  real,parameter :: dt_SI = omega / omega_SI, c_u_SI = dx_SI / dt_SI, c_s_SI = c_u_SI /3**(1/2.)
end module declare

subroutine make_Lists(PL_Liqx,PL_Liqy,&
                      LL_Liqx,LL_Liqy,LL_Liqi, &
                      LL_Borx,LL_Bory,LL_Bori, &
                      LL_totx,LL_toty,LL_toti, &
                      VelBorx,VelBory, &
                      Li2Source,Li2ibar,Po2Li, &
                      Bot2Bor,x02Po, &
                      nPL_Liq,nLL_Liq,nLL_Bor,nLL_tot,nBot2Bor)
  use declare
  implicit none
  integer :: PL_Liqx(1:nPo),PL_Liqy(1:nPo) 
  integer :: LL_Liqx(1:nPoLi),LL_Liqy(1:nPoLi),LL_Liqi(1:nPoLi)
  integer :: LL_Borx(1:nPoLi),LL_Bory(1:nPoLi),LL_Bori(1:nPoLi)
  integer :: LL_totx(1:nPoLi),LL_toty(1:nPoLi),LL_toti(1:nPoLi)
  complex :: VelBorx(1:nPoLi),VelBory(1:nPoLi)
  integer :: Li2Source(1:nPoLi),Li2ibar(1:nPoLi),Po2Li(1:nPo,0:ndir)
  integer :: Bot2Bor(1:nx),x02Po(1:ny)
  integer :: nPL_Liq,nLL_Liq,nLL_Bor,nLL_tot,nBot2Bor
  integer :: x,y,xmd,ymd,i,il,ip,il2,x2,y2,i2

  nPL_Liq = 0; nLL_Liq = 0; nLL_Bor = 0; nBot2Bor = 0
  DO y = 1,ny
    x02Po(y) = nPL_Liq + 1
    DO x = 1,nx
      IF ((y.le.ny).and.(y.ge.1).and.((x-x_sphere)**2 +(y-y_sphere)**2.gt.rad_sphere**2)) THEN
        nPL_Liq = nPL_Liq + 1
        PL_Liqx(nPL_Liq) = x
        PL_Liqy(nPL_Liq) = y
        DO i = 0,ndir
          xmd = MOD(nx+x-dirx(i)-1,nx)+1
          ymd = y-diry(i)
          IF ((ymd.le.ny).and.(ymd.ge.1).and.((xmd-x_sphere)**2 +(ymd-y_sphere)**2.gt.rad_sphere**2)) THEN
            nLL_Liq = nLL_Liq + 1
            LL_Liqx(nLL_Liq) = x
            LL_Liqy(nLL_Liq) = y
            LL_Liqi(nLL_Liq) = i
          ELSE
            nLL_Bor = nLL_Bor + 1
            LL_Borx(nLL_Bor) = x
            LL_Bory(nLL_Bor) = y
            LL_Bori(nLL_Bor) = i
            IF(ymd.gt.ny) THEN
              VelBorx(nLL_Bor) = complex(0,0)
              VelBory(nLL_Bor) = complex(0,0)
            ELSEIF(((xmd-x_sphere)**2 +(ymd-y_sphere)**2.lt.rad_sphere**2).OR.(ymd.lt.1)) THEN
              nBot2Bor = nBot2Bor + 1
              Bot2Bor(nBot2Bor) = nLL_Bor
              VelBorx(nLL_Bor) = complex(uxWall,0)
              VelBory(nLL_Bor) = complex(0     ,0)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDDO
  nLL_tot = nLL_Liq + nLL_Bor
  DO il = 1,nLL_Liq
    LL_totx(il) = LL_Liqx(il)
    LL_toty(il) = LL_Liqy(il)
    LL_toti(il) = LL_Liqi(il)
  ENDDO
  DO il = 1,nLL_Bor
    LL_totx(il+nLL_Liq) = LL_Borx(il)
    LL_toty(il+nLL_Liq) = LL_Bory(il)
    LL_toti(il+nLL_Liq) = LL_Bori(il)
  ENDDO

  DO il = 1,nLL_Liq
    x = LL_Liqx(il)
    y = LL_Liqy(il)
    i = LL_Liqi(il)
    xmd = MOD(nx+x-dirx(i)-1,nx)+1
    ymd = y-diry(i)
    DO il2 = 1,nLL_tot
      x2 = LL_totx(il2)
      y2 = LL_toty(il2)
      i2 = LL_toti(il2)
      IF ((xmd == x2).and.(ymd == y2).and.(i == i2)) THEN
        Li2Source(il) = il2
      ENDIF
    ENDDO
  ENDDO

  DO il = 1,nLL_Bor
    x = LL_Borx(il)
    y = LL_Bory(il)
    i = LL_Bori(il)
    DO il2 = 1,nLL_tot
      x2 = LL_totx(il2)
      y2 = LL_toty(il2)
      i2 = LL_toti(il2)
      IF ((x == x2).and.(y == y2).and.(i2 == no_slip(i))) THEN
        Li2ibar(il) = il2
      ENDIF
    ENDDO
  ENDDO

  print *,'nPL_Liq',nPL_Liq
  DO ip = 1,nPL_Liq
    x = PL_Liqx(ip)
    y = PL_Liqy(ip)
    DO i = 0,ndir
      DO il = 1,nLL_tot
        x2 = LL_totx(il)
        y2 = LL_toty(il)
        i2 = LL_toti(il)
        IF ((x == x2).and.(y == y2).and.(i2 == i)) THEN
          Po2Li(ip,i) = il
        ENDIF
      ENDDO
    ENDDO
  ENDDO
END subroutine make_Lists

program repeat
  use declare
  implicit none
  complex :: dfc
  complex :: ux(1:nPo),uy(1:nPo),dr(1:nPo)
  complex :: h1(1:nPoLi),h2(1:nPoLi)
  complex :: drdis(1:nx,1:ny),uxdis(1:nx,1:ny),uydis(1:nx,1:ny)
  integer :: PL_Liqx(1:nPo),PL_Liqy(1:nPo)
  integer :: LL_Liqx(1:nPoLi),LL_Liqy(1:nPoLi),LL_Liqi(1:nPoLi)
  integer :: LL_Borx(1:nPoLi),LL_Bory(1:nPoLi),LL_Bori(1:nPoLi)
  integer :: LL_totx(1:nPoLi),LL_toty(1:nPoLi),LL_toti(1:nPoLi)
  complex :: VelBorx(1:nPoLi),VelBory(1:nPoLi),ForceBorx(1:nPoLi),ForceBory(1:nPoLi)
  integer :: Li2Source(1:nPoLi),Li2ibar(1:nPoLi),Po2Li(1:nPo,0:ndir)
  integer :: Bot2Bor(1:nx),x02Po(1:ny)
  integer :: nPL_Liq,nLL_Liq,nLL_Bor,nBot2Bor,nLL_tot
  integer :: x,y,xmd,ymd,i,ibar,il,ip,ilB,il2,x2,y2,i2

  DO ip = 1,nPo
    dr(ip) = complex(0,0)
    ux(ip) = complex(0,0)
    uy(ip) = complex(0,0)
  ENDDO
!  open(unit = 20,file='complex_freq.asc1',action='write')
  ! write (20,*) NSTEPS

  call make_Lists(PL_Liqx,PL_Liqy,&
                      LL_Liqx,LL_Liqy,LL_Liqi, &
                      LL_Borx,LL_Bory,LL_Bori, &
                      LL_totx,LL_toty,LL_toti, &
                      VelBorx,VelBory, &
                      Li2Source,Li2ibar,Po2Li, &
                      Bot2Bor,x02Po, &
                      nPL_Liq,nLL_Liq,nLL_Bor,nLL_tot,nBot2Bor)

  call init_equilibrium(h1,dr,ux,uy,Po2Li)
  do i = 1,NSTEPS
    call stream(h1,h2,VelBorx,VelBory,ForceBorx,ForceBory,nLL_Liq,nLL_Bor,Li2Source,LL_Bori,LL_toti,Li2ibar)
    call compute_dr_u(h2,dr,ux,uy   ,nPL_Liq,Po2Li)
    call collide(     h2,dr,ux,uy,h1,nPL_Liq,Po2Li)
    call calc_dfc(ForceBorx,nBot2Bor,dfc)
    print *, NSTEPS,real(dfc),aimag(dfc)
!    write (20,*) real(dfc),aimag(dfc)
  end do
  ! close(unit = 20)

  ! open(unit = 30,file='uxuy.asc1',action='write')
  ! write(30,*)nx,ny
  ! DO x=1,nx
  !   DO y=1,ny
  !     write(30,*)x,y,uxdis(x,y),uydis(x,y)
  !   ENDDO
  ! ENDDO
  ! close(unit = 30)
  print *,'we are done'

end program repeat

subroutine init_equilibrium(h1,dr,ux,uy,Po2Li)
  use declare
  implicit none
  complex :: cidotu
  integer :: ip,i,il
  complex :: h1(1:nPoLi)
  complex :: dr(1:nPo),ux(1:nPo),uy(1:nPo)
  integer :: Po2Li(1:nPo,0:ndir)
  DO ip = 1,nPo
    DO i = 0,ndir
      il = Po2Li(ip,i)
      cidotu = dirx(i)*ux(ip)+diry(i)*uy(ip)
      h1(il) = wi(i)*dr(ip)*(1+3*cidotu)
    END DO
  END DO
END subroutine init_equilibrium

subroutine stream(h_src,h_dst,VelBorx,VelBory,ForceBorx,ForceBory,nLL_Liq,nLL_Bor,Li2Source,LL_Bori,LL_toti,Li2ibar)
  use declare
  implicit none
  complex :: h_src(nPoLi),h_dst(nPoLi)
  complex :: VelBorx(1:nPoLi),VelBory(1:nPoLi),ForceBorx(1:nPoLi),ForceBory(1:nPoLi)
  integer :: il,ilB,i,ibar,nLL_Liq,nLL_Bor
  integer :: Li2Source(1:nPoLi),Li2ibar(1:nPoLi)
  integer :: LL_Bori(1:nPoLi),LL_toti(1:nPoLi)
  complex :: cidotu

  DO il = 1,nLL_Liq
    h_dst(il) = h_src(Li2Source(il))
  ENDDO

  DO ilB = 1,nLL_Bor
    il = ilB + nLL_Liq
    i = LL_Bori(il)
    ibar = LL_toti(Li2ibar(ilB))
    cidotu = dirx(ibar)*VelBorx(ilB) + diry(ibar)*VelBory(ilB)
    h_dst(il) = h_src(Li2ibar(ilB))- 6 * wi(ibar) * cidotu      !# Eq. 5.26
    ForceBorx(ilB) = (h_dst(il) + h_dst(Li2ibar(ilB)))*dirx(i)*wi(i)    !# Eqs. 5.79, 5.80
    ForceBory(ilB) = (h_dst(il) + h_dst(Li2ibar(ilB)))*diry(i)*wi(i)    !# Eqs. 5.79, 5.80
  ENDDO
END subroutine stream

subroutine compute_dr_u(h_in,dr,ux,uy,nPL_Liq,Po2Li)
  use declare
  implicit none
  complex :: h_in(1:nPoLi)
  complex :: ux(1:nPo),uy(1:nPo),dr(1:nPo)
  integer :: nPL_Liq,ip,il,i
  integer :: Po2Li(1:nPo,0:ndir)
  DO ip = 1,nPL_Liq
    dr(ip) = complex(0,0)
    ux(ip) = complex(0,0)
    uy(ip) = complex(0,0)
    DO i = 0,ndir
      il = Po2Li(ip,i)
      dr(ip) = dr(ip) + h_in(il)
      ux(ip) = ux(ip) + h_in(il)*dirx(i)
      uy(ip) = uy(ip) + h_in(il)*diry(i)
    ENDDO
  ENDDO
END subroutine compute_dr_u

subroutine collide(h_in,dr,ux,uy,h_out,nPL_Liq,Po2Li)
  use declare
  implicit none
  complex :: h_in(1:nPoLi),h_out(1:nPoLi),cidotu,heq,ux(1:nPo),uy(1:nPo),dr(1:nPo)
  integer :: nPL_Liq,ip,i,il,Po2Li(1:nPo,0:ndir)
  DO ip = 1,nPL_Liq
    DO i = 0,ndir
      cidotu = dirx(i)*ux(ip) + diry(i)*uy(ip)
      heq = dr(ip) + 3.*cidotu
      il = Po2Li(ip,i)
      h_out(il) = omtauinv * h_in(il) + tauinv * heq + complex(0.,1.)*omega * h_in(il) - omega**2./2. * heq   ! # eq. 3.3
    ENDDO
  ENDDO
END subroutine collide

subroutine calc_dfc(ForceBorx,nBot2Bor,dfc)
  use declare
  implicit none
  complex :: TanForce,TanStress,ZL,Zq,dfc,ForceBorx(1:nPoLi)
  integer :: il,nBot2Bor

  TanForce = complex(0.,0.)
  DO il = 1,nBot2Bor
    TanForce = TanForce + ForceBorx(il)
    TanStress = TanForce / nx
  ENDDO
  ZL = TanStress/uxWall *c_U_SI/ delta
  Zq = (rho0/rho_Liq_SI)*(nu/nu_SI)*(delta_SI/delta)*Zq_SI
  dfc = f0_SI* complex(0.,1.)/ Pi * ZL/Zq
  dfc = -CONJG(dfc)
end subroutine calc_dfc
