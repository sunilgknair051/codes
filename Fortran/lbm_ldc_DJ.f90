program lbm
! -------------------------------------------------------------
! scheme for the solution of an incompressible flow.
! D2Q9 velocity LBM-BGK model notation of velocity stencil:
!     6   2   5
!       \ | /
!     3 - 0 - 1
!       / | \
!     7   4   8
! boundary conditions are set for "lid driven cavity" problem classical approach, no performance optimization
! compile with -fast (on IA64  -ipo, -O3, -no-prec-div, -static, and -xP )

  implicit none

  integer it,itmax,imax,jmax,nspeed,itstep
  integer i,j,n,ie,iw,jn,js,nv
!  parameter(imax=100,jmax=100,nspeed=8)
!  real f(0:nspeed,1:imax,1:jmax),fn(0:nspeed,1:imax,1:jmax)
!  real fe(0:nspeed),u,v,r,u2,v2,upv,umv,upv2,umv2,t0,t1,t2
  real(kind=8),allocatable ::  f(:,:,:),fn(:,:,:),fe(:)
  !integer, allocatable, dimension (:) :: i1
  real u,v,r,u2,v2,upv,umv,upv2,umv2,t0,t1,t2
  real vis,omega,eomega,usq,dens,lidvel,t1rl,t2rl,t1x,t2x,force
  real time0,time1

!     input parameter from stdin

  !open(66,file='mlups_nonopt.dat',form='formatted')
  !open(15,file='run.dat',form='formatted')
  !read(15,*) vis
  !read(15,*) force
  !read(15,*) nmax
  !allocate (i1(1:nmax))
  !read(15,*) (i1(n),n=1,nmax)
  !close(15)

  nspeed = 8
  omega  = 1./(0.5 + 3.0*vis )
  eomega = 1.0-omega
  dens   = 1.0
  t1x    = force/6.0
  vis = 1/6
  force = 0.00001

!     init with equilibrium distribution function for zero velocity

  t0  = dens*4.d0/9.d0
  t1  = dens*1.d0/9.d0
  t2  = dens*1.d0/36.d0

  ! do nv=1,nmax
  !   imax = i1(nv)
  !   jmax = i1(nv)
     imax = 10
     jmax = 10
     allocate (f(0:nspeed,1:imax,1:jmax),fn(0:nspeed,1:imax,1:jmax),fe(0:nspeed))
     do j=1,jmax
        do i=1,imax
           f(0,i,j)=t0
           f(1,i,j)=t1
           f(2,i,j)=t1
           f(3,i,j)=t1
           f(4,i,j)=t1
           f(5,i,j)=t2
           f(6,i,j)=t2
           f(7,i,j)=t2
           f(8,i,j)=t2
        enddo
     enddo
!     start main iteration loop
     itmax = 4e6*30/(imax**2)
     itmax = 10000
     itstep = itmax/10
     call cpu_time(time0)
     do it=1,itmax
        if(mod(it,itstep).eq.0)print*,it
!     propagation step, assuming periodic boundary conditions
        do i=1,imax
           do j=1,jmax
              ie = mod(i,imax) + 1
              iw = imax - mod(imax+1-i,imax)
              jn = mod(j,jmax) + 1
              js = jmax - mod(jmax+1-j,jmax)
              fn(1,ie ,j ) = f(1,i,j)
              fn(2,i  ,jn) = f(2,i,j)
              fn(3,iw ,j ) = f(3,i,j)
              fn(4,i  ,js) = f(4,i,j)
              fn(5,ie ,jn) = f(5,i,j)
              fn(6,iw ,jn) = f(6,i,j)
              fn(7,iw ,js) = f(7,i,j)
              fn(8,ie ,js) = f(8,i,j)
              fn(0,i  ,j ) = f(0,i,j)
           enddo
        enddo
!     boundary conditions: bounce back
!     lower wall j=1
        j=1
        do i=2,imax-1
           f(0,i,j) =  fn(0,i,j)
           f(1,i,j) =  fn(3,i,j)
           f(2,i,j) =  fn(4,i,j)
           f(3,i,j) =  fn(1,i,j)
           f(4,i,j) =  fn(2,i,j)
           f(5,i,j) =  fn(7,i,j)
           f(6,i,j) =  fn(8,i,j)
           f(7,i,j) =  fn(5,i,j)
           f(8,i,j) =  fn(6,i,j)
        enddo
!     upper wall j=jmax
        j=jmax
        do i=2,imax-1
           f(0,i,j) =  fn(0,i,j)
           f(1,i,j) =  fn(3,i,j)
           f(2,i,j) =  fn(4,i,j)
           f(3,i,j) =  fn(1,i,j)
           f(4,i,j) =  fn(2,i,j)
           f(5,i,j) =  fn(7,i,j)
           f(6,i,j) =  fn(8,i,j)
           f(7,i,j) =  fn(5,i,j)
           f(8,i,j) =  fn(6,i,j)
        enddo
!     west wall i=1
        i=1
        do j=2,jmax-1
           f(0,i,j) =  fn(0,i,j)
           f(1,i,j) =  fn(3,i,j)
           f(2,i,j) =  fn(4,i,j)
           f(3,i,j) =  fn(1,i,j)
           f(4,i,j) =  fn(2,i,j)
           f(5,i,j) =  fn(7,i,j)
           f(6,i,j) =  fn(8,i,j)
           f(7,i,j) =  fn(5,i,j)
           f(8,i,j) =  fn(6,i,j)
        enddo
!     east wall: i=imax
        i=imax
        do j=2,jmax-1
           f(0,i,j) =  fn(0,i,j)
           f(1,i,j) =  fn(3,i,j)
           f(2,i,j) =  fn(4,i,j)
           f(3,i,j) =  fn(1,i,j)
           f(4,i,j) =  fn(2,i,j)
           f(5,i,j) =  fn(7,i,j)
           f(6,i,j) =  fn(8,i,j)
           f(7,i,j) =  fn(5,i,j)
           f(8,i,j) =  fn(6,i,j)
        enddo
!     north wall: add source term to drive flow in x-direction
        j=jmax
        do i=2,imax-1
           f(7,i,j) =  f(7,i,j) -  t1x
           f(8,i,j) =  f(8,i,j) +  t1x
        enddo
!     relaxation step
        t0  = 4.d0/9.d0
        t1  = 1.d0/9.d0
        t2  = 1.d0/36.d0
        do j = 2,jmax-1
           do i = 2,imax-1
              u = fn(1,i,j) + fn(5,i,j) + fn(8,i,j)  - fn(6,i,j) - fn(3,i,j) - fn(7,i,j)
              v = fn(5,i,j) + fn(2,i,j) + fn(6,i,j)  - fn(7,i,j) - fn(4,i,j) - fn(8,i,j)
              r = fn(0,i,j) + fn(1,i,j) + fn(2,i,j)  + fn(3,i,j) + fn(4,i,j) + fn(5,i,j) + fn(6,i,j) + fn(7,i,j) + fn(8,i,j)

              u = u/r
              v = v/r
              t1rl = t1*r
              t2rl = t2*r
              u2  = u*u
              v2  = v*v
              usq = 1.5d0*(u2+v2)

              upv = 3.d0*(u+v)
              umv = 3.d0*(u-v)
              upv2 = 0.5d0*upv*upv
              umv2 = 0.5d0*umv*umv

              fe(0) = t0*r* (1.d0 - usq )
              fe(1) = t1rl * (1.d0 + 3.d0*u + 4.5d0*u2  - usq )
              fe(2) = t1rl * (1.d0 + 3.d0*v + 4.5d0*v2  - usq )
              fe(3) = t1rl * (1.d0 - 3.d0*u + 4.5d0*u2  - usq )
              fe(4) = t1rl * (1.d0 - 3.d0*v + 4.5d0*v2  - usq )

              fe(5) = t2rl * (1.d0 + upv + upv2 - usq )
              fe(6) = t2rl * (1.d0 - umv + umv2 - usq )
              fe(7) = t2rl * (1.d0 - upv + upv2 - usq )
              fe(8) = t2rl * (1.d0 + umv + umv2 - usq )

              do n=0,nspeed
                 f(n,i,j) = eomega*fn(n,i,j) + omega*fe(n)
              enddo

           enddo
        enddo

     enddo
! --------------------------------------------------------------
!     end main interation loop
! --------------------------------------------------------------

     call cpu_time(time1)
     write(6,*)
     write(6,*) ' Run # : '
     write(6,*) ' CPU Time used :',time1-time0
     write(6,*) ' Iterations    :',itmax
     write(6,*) ' Grid          :',imax,jmax
     write(6,*) ' MLUP/s        :',1e-6*imax*jmax*itmax/(time1-time0)

     write(66,*) imax,1e-6*imax*jmax*itmax/(time1-time0)

!     goto 1111

! --------------------------------------------------------------
!     output results in tecplot format
! --------------------------------------------------------------

     open(11,file='tecout.tec',form='formatted')
     write(11,*) 'variables = "x", "y", "u" , "v" , "r" '
     write(11,*) 'zone i=',imax,' j=',jmax,' f=point'
     write(11,*) ''
     do j=1,jmax
        do i=1,imax
           r =  f(0,i,j) + f(1,i,j) + f(2,i,j) + f(3,i,j) + f(4,i,j) + f(5,i,j) + f(6,i,j) + f(7,i,j) + f(8,i,j)
           u =  f(1,i,j) + f(5,i,j) + f(8,i,j) - f(6,i,j) - f(3,i,j) - f(7,i,j)
           v =  f(5,i,j) + f(2,i,j) + f(6,i,j) - f(7,i,j) - f(4,i,j) - f(8,i,j)
           if(i.eq.1) then
              u=0.
              v=0.
              r=1.
           endif
           if(i.eq.imax)then
              u=0.
              v=0.
              r=1.
           endif
           if(j.eq.1) then
              u=0.
              v=0.
              r=1.
           endif
           if(j.eq.jmax) then
              u=force
              v=0.
              r=1.
           endif
           write(11,1) float(i),float(j),u/r,v/r,r
        enddo
     enddo
     close(11)

1111 deallocate (f,fn,fe)

! enddo

  !deallocate (i1)
  close(66)
! --------------------------------------------------------------
!     end program lbm
! --------------------------------------------------------------

  stop
1 format (5(1x,1pe12.5))
end program lbm
