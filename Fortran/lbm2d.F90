! !     ------------------------------------------------------------------
! !       LATTICE BOLTZMANN METHOD FOR SIMULATING DRIVEN CAVITY FLOW
! !    ------------------------------------------------------------------
! !     * 2-Dimensional geometry - Uniform grid along both X and Z
! !     * Bhatnagar-Gross-Krook (BGK) FINITE DIFFERENCE discretization
! !     * 9-velocity model { 0 1 2 3 4 5 6 7 8 }
! !     * Unsteady, incompressible (almost !), viscous flow
! !     * Bounce-back boundary conditions on stationary walls
! !     * Postprocessing using MATLAB
! !
! !   BENCHMARK PROBLEMS
! !
! !       1. The 2-Dimensional driven cavity
! !
! !    Written By    : Abhijit S. Joshi
! !     Last modified : Thu 25 Oct 2001 at 12:30 hrs
     ! Remarks       : Studying the effect of Mach number on accuracy
     ! References : [1] Simulation of Cavity Flow by the Lattice
                     ! Boltzmann method - Hou et.al.(1995)
                      ! Journal of Computational Physics,118,329-347
                      ! [2] Accuracy and efficiency study of LBM for
                     ! steady state flow simulations - Lai et.al.(2001)
                      ! Numerical Heat Transfer(B),39,21-43
     ! ------------------------------------------------------------------
                             ! Main Program
     ! ------------------------------------------------------------------
      INCLUDE 'gridsize.inc'         ! contains { nx , nz }
      INCLUDE 'common.inc'           ! contains declarations & common

      Re = 100.0D0                   ! Reynolds number

      CALL grid                      ! set up grid for the problem
      CALL initial                   ! specify initial conditions
      CALL time_integration          ! Evolve the "f" distributions
      CALL write_results             ! create output files

      END
! C     ------------------------------------------------------------------
! C                       * SETTING UP THE GRID *
! C
! C      nx = no. of nodes along the X direction (1,2,3,....nx)
! C      nz = no. of nodes along the Z direction (1,2,3,....nz)
! C      ( x(i), z(k) )  =  X-Z coordinates of node (i,k)
! C
! C     ------------------------------------------------------------------
      SUBROUTINE grid
      INCLUDE 'gridsize.inc'
      INCLUDE 'common.inc'

      WRITE(6,*)
      WRITE(6,*)'Grid generation........'
      WRITE(6,*)

! c     Calculating grid spacing along X and Z.

      xl = 1.0D0             ! domain length along X in lattice units
      zl = 1.0D0             ! domain length along Z in lattice units

      dx = xl/FLOAT(nx-1)    ! lattice spacing along X
      dz = zl/FLOAT(nz-1)    ! lattice spacing along Z
      dt = dx                ! time step

! c     Calculating the X and Z coordinates of the node points

      x(1) = 0.0D0
      x(nx) = xl
      DO i=2,nx-1
       x(i) = x(i-1) + dx
      ENDDO

      z(1) = 0.0D0
      z(nz) = zl
      DO k=2,nz-1
       z(k) = z(k-1) + dz
      ENDDO

! c     calculating the "direction vectors" for LBM. The 9-velocity model
! c     is used here with the numbers running from {0,1,2,3,4,5,6,7,8}

      ex(0) =  0.0D0
      ez(0) =  0.0D0

      ex(1) =  1.0D0
      ez(1) =  0.0D0

      ex(2) =  1.0D0
      ez(2) =  1.0D0

      ex(3) =  0.0D0
      ez(3) =  1.0D0

      ex(4) = -1.0D0
      ez(4) =  1.0D0

      ex(5) = -1.0D0
      ez(5) =  0.0D0

      ex(6) = -1.0D0
      ez(6) = -1.0D0

      ex(7) =  0.0D0
      ez(7) = -1.0D0

      ex(8) =  1.0D0
      ez(8) = -1.0D0

      END
! C     ------------------------------------------------------------------
! C     This is the initialization for the 2-dimensional Lid driven cavity
! C     flow. The velocity is set to zero everywhere inside the cavity and
! C     density is set equal to unity. Later, velocity at the lid is set
! C     to UINF, which is designed to match the desired Reynolds number.
! C       As far as the LBM is concerned, the equilibrium distribution
! C     functions are calculated here based on the initial known velocity
! C     field both on the boundaries and interior nodes.
! C     ------------------------------------------------------------------
      SUBROUTINE initial
      INCLUDE 'gridsize.inc'
      INCLUDE 'common.inc'

      WRITE(6,*)
      WRITE(6,*)'Initialization........'
      WRITE(6,*)

! c     Initially the fluid inside the cavity is at rest, and at t = 0
! c     the top layer starts moving at speed "U" to the right. The values
! c     of the properties are adjusted to obtain the desired Reynolds no.

      densf  = 1.0D0                           ! reference density
      tau = 0.60D0                             ! relaxation time
      viskin = (2.0D0*tau - 1.0D0)*dx**2/(6.0D0*dt)       ! kinematic viscosity
      UINF = viskin * Re / xl                  ! suitable lid velocity
      cs = dx/(dt*SQRT(3.0D0))                 ! speed of sound

      WRITE(6,*)'Parameter values used in the simulation:'
      WRITE(6,*)
      WRITE(6,*)'Reynolds number (Re) = ',Re
      WRITE(6,*)'Domain size (xl zl)  = ',xl,zl
      WRITE(6,*)'Grid size   (nx nz)  = ',nx,nz
      WRITE(6,*)'grid spacing dx = dz = ',dx
      WRITE(6,*)'Time step (dt)       = ',dt
      WRITE(6,*)'Reference density    = ',densf
      WRITE(6,*)'Kinematic viscosity  = ',viskin
      WRITE(6,*)'Relaxation time      = ',tau
      WRITE(6,*)'Lid velocity (UINF)  = ',UINF
      WRITE(6,*)'Speed of sound (cs)  = ',cs
      WRITE(6,*)'Mach number          = ',UINF/cs
      WRITE(6,*)

      DO k=1,nz
      DO i=1,nx
       den(i,k) = densf           ! set fluid density
       u(i,k) = 0.0D0             ! set fluid velocity = 0
       w(i,k) = 0.0D0
      ENDDO
      ENDDO

      DO i=1,nx                   ! top lid moves to the right
       u(i,nz) = UINF
      ENDDO

! c     weight functions used in the equilibrium distribution

      alpha(0) = 4.0D0/9.0D0
      DO id = 1,7,2
       alpha(id) = 1.0D0/9.0D0
      ENDDO
      DO id = 2,8,2
       alpha(id) = 1.0D0/36.0D0
      ENDDO

! c     The equilibrium particle distributions are computed based on
! c     the initial velocity distributions.

      CALL equilibrium_distribution

! c     initial equilibrium distributions are stored here

      DO id=0,8
      DO i=1,nx
      DO k=1,nz
       feq_init(i,k,id) = feq(i,k,id)
      ENDDO
      ENDDO
      ENDDO

! c     particle distributions are assigned equilibrium values

      DO id=0,8
      DO k=1,nz
      DO i=1,nx
       f(i,k,id)     = feq(i,k,id)
       fprev(i,k,id) = f(i,k,id)
      ENDDO
      ENDDO
      ENDDO

! c     Does the equilibrium distribution give back the initial density
! c     and velocity everywhere in the domain ?

      WRITE(6,*)'Checking density & velocity after initialization'

      CALL density          ! recompute the density    den(i,k)
      CALL velocity         ! recompute the velocity   u(i,k), w(i,k)

      d_den = 0.0D0
      d_vex = 0.0D0
      d_vez = 0.0D0

      DO i=1,nx
       d_den = MAX( DABS(den(i,nz) - densf), d_den)
       d_vex = MAX( DABS(u(i,nz) - UINF), d_vex)
       d_vez = MAX( DABS(w(i,nz) - 0.0D0), d_vez)
      ENDDO

      DO i=1,nx
      DO k=1,nz-1
       d_den = MAX( DABS(den(i,k) - densf), d_den)
       d_vex = MAX( DABS(u(i,k) - 0.0D0), d_vex)
       d_vez = MAX( DABS(w(i,k) - 0.0D0), d_vez)
      ENDDO
      ENDDO

      WRITE(6,*)
      WRITE(6,*)'Max.change in density    = ',d_den
      WRITE(6,*)'              x-velocity = ',d_vex
      WRITE(6,*)'              z-velocity = ',d_vez
      WRITE(6,*)

      END
! C     ------------------------------------------------------------------
! C     Calculate the Equilibrium particle distribution functions at all
! C     node points based on the current velocity & density fields.
! C     ------------------------------------------------------------------
      SUBROUTINE equilibrium_distribution
      INCLUDE 'gridsize.inc'
      INCLUDE 'common.inc'

      DO id = 0,8
      DO k=1,nz
      DO i=1,nx
       feq(i,k,id) = den(i,k) * alpha(id) * ( 1.0 &
                  + 3.0 * (ex(id)*u(i,k) + ez(id)*w(i,k)) &
                  + 4.5 * (ex(id)*u(i,k) + ez(id)*w(i,k))**2 &
                  - 1.5 * (u(i,k)**2 + w(i,k)**2) )
      ENDDO
      ENDDO
      ENDDO

      END
! C     ------------------------------------------------------------------
! C                        TIME INTEGRATION ROUTINE
! C
! C       The "f" distributions {0,1,2,3,4,5,6,7,8} are advanced in time
! C     based on the discrete LGBK model. (Explicit method; No iteration)
! C     ------------------------------------------------------------------
      SUBROUTINE time_integration
      INCLUDE 'gridsize.inc'
      INCLUDE 'common.inc'

      time = 0.0D0

      WRITE(6,*)
      WRITE(6,*)'TIME INTEGRATION BEGINS'
      WRITE(6,*)

! c     epsnss = 1.0D-6             ! criterion for steady state

      DO itime=1,100000           ! time-integration begins

       time = time + dt           ! advance time to the next level

       CALL collision             ! evaluate the R.H.S.
       CALL streaming             ! update the node points
       CALL boundary_conditions   ! bounce back scheme implemented
       CALL density               ! construct the density field
       CALL velocity              ! construct the velocity field

       CALL check_steady_state(dfdt)
       IF(MOD(itime,100).EQ.0) WRITE(6,100)INT(itime),time,dfdt
       !IF(dfdt.LT.epsnss) GOTO 200

      ENDDO                ! time-integration ends

      WRITE(6,*)'WARNING:- Terminating before steady state'
      RETURN

      WRITE(6,*)'Steady state reached. Terminating computation'
      WRITE(6,*)

      FORMAT(1X,'time step = ',I8,3X,'t = ',E12.5,3X,'dfdt = ',E12.5)
      END
! C     ------------------------------------------------------------------
! C           Calculate the R.H.S. of the discrete LBGK equation
! C     ------------------------------------------------------------------
      SUBROUTINE collision
      INCLUDE 'gridsize.inc'
      INCLUDE 'common.inc'

      CALL equilibrium_distribution

      DO id = 0,8
      DO k=1,nz
      DO i=1,nx
       omega(i,k,id) = - ( fprev(i,k,id) - feq(i,k,id) ) / tau
      ENDDO
      ENDDO
      ENDDO

      END
! C     ------------------------------------------------------------------
! C                Update the "f" fields at all node points
! C     ------------------------------------------------------------------
      SUBROUTINE streaming
      INCLUDE 'gridsize.inc'
      INCLUDE 'common.inc'

      DO id=0,8
      DO k=1,nz
      DO i=1,nx
       idelx = INT( ex(id) )
       idelz = INT( ez(id) )
       f(i+idelx,k+idelz,id) = fprev(i,k,id) + omega(i,k,id)
      ENDDO
      ENDDO
      ENDDO

      END
! C     ------------------------------------------------------------------
! C                          Bounce Back Scheme
! C     ------------------------------------------------------------------
      SUBROUTINE boundary_conditions
      INCLUDE 'gridsize.inc'
      INCLUDE 'common.inc'

! c     top distribution functions set to their initial equilibrium values

      DO id=0,8
      DO i=1,nx
       f(i,nz,id) = feq_init(i,nz,id)
      ENDDO
      ENDDO

! c     Bounce back scheme is used to determine the distribution functions
! c     on walls with the "no-slip" velocity boundary conditions

! c     Left wall ( f8 f1 f2 )

      DO k=2,nz-1
       f(1,k,8) = f(1,k,4)
       f(1,k,1) = f(1,k,5)
       f(1,k,2) = f(1,k,6)
      ENDDO

! c     Right wall ( f4 f5 f6 )

      DO k=2,nz-1
       f(nx,k,4) = f(nx,k,8)
       f(nx,k,5) = f(nx,k,1)
       f(nx,k,6) = f(nx,k,2)
      ENDDO

! c     Bottom wall ( f2 f3 f4 )

      DO i=2,nx-1
       f(i,1,2) = f(i,2,6)
       f(i,1,3) = f(i,2,7)
       f(i,1,4) = f(i,2,8)
      ENDDO

! c     Bottom left corner ( f8 f1 f2 f3 f4 )

      f(1,1,1) = f(1,1,5)
      f(1,1,2) = f(1,1,6)
      f(1,1,3) = f(1,1,7)
      f(1,1,4) = ( den(1,1) - (f(1,1,0) + f(1,1,1) + f(1,1,2) + f(1,1,3) &
                           +  f(1,1,5) + f(1,1,6) + f(1,1,7)) ) * 0.5
      f(1,1,8) = f(1,1,4)

! c     Bottom right corner ( f2 f3 f4 f5 f6 )

      f(nx,1,3) = f(nx,1,7)
      f(nx,1,4) = f(nx,1,8)
      f(nx,1,5) = f(nx,1,1)
      f(nx,1,2) = ( den(nx,1) - (f(nx,1,0) + f(nx,1,1) + f(nx,1,3) &
           +  f(nx,1,4) + f(nx,1,5) + f(nx,1,7)+  f(nx,1,8)) ) * 0.5
      f(nx,1,6) = f(nx,1,2)

      END
! C     ------------------------------------------------------------------
! C                Compute the density at all node points
! C     ------------------------------------------------------------------
      SUBROUTINE density
      INCLUDE 'gridsize.inc'
      INCLUDE 'common.inc'

      DO k=1,nz
      DO i=1,nx
       sum = 0.0D0
       DO id = 0,8
        sum = sum + f(i,k,id)
       ENDDO
       den(i,k) = sum
      ENDDO
      ENDDO

      END
! C     ------------------------------------------------------------------
! C                Compute the velocity at all node points
! C     ------------------------------------------------------------------
      SUBROUTINE velocity
      INCLUDE 'gridsize.inc'
      INCLUDE 'common.inc'

      DO k=1,nz
      DO i=1,nx
       sumx = 0.0D0
       sumz = 0.0D0
       DO id = 0,8
        sumx = sumx + f(i,k,id)*ex(id)
        sumz = sumz + f(i,k,id)*ez(id)
       ENDDO
       u(i,k) = sumx / den(i,k)
       w(i,k) = sumz / den(i,k)
      ENDDO
      ENDDO

      END
! C     ------------------------------------------------------------------
! C                Compute the pressure at all node points
! C     ------------------------------------------------------------------
      SUBROUTINE pressure
      INCLUDE 'gridsize.inc'
      INCLUDE 'common.inc'

! c     In the LBM, the density is not really a constant even for
! c     incompressible flows, and the pressure is proportional to the
! c     density. But in this problem, the pressure is not required in the
! c     calculation.

      C_s = 1.0D0/DSQRT(3.0D0)       ! speed of sound in the fluid
      DO i=1,nx
      DO k=1,nz
       p(i,k) = C_s**2 * den(i,k)
      ENDDO
      ENDDO

      END
! C     ------------------------------------------------------------------
! C             Checking if the system has reached steady state
! C     ------------------------------------------------------------------
      SUBROUTINE check_steady_state(dfdt)
      INCLUDE 'gridsize.inc'
      INCLUDE 'common.inc'

      dfdt = 0.0D0
      DO k=1,nz
      DO i=1,nx
      DO id = 0,8
       chf = DABS( f(i,k,id) - fprev(i,k,id) )
       dfdt = MAX( dfdt, chf )
       fprev(i,k,id) = f(i,k,id)
      ENDDO
      ENDDO
      ENDDO

      END
