program heat

   implicit none
   include 'mpif.h'

   ! Sizes for discretization
   integer size_x, size_y, size_z, me, x_domains, y_domains, z_domains
   integer size_x_glo, size_y_glo, size_z_glo

   ! Arrays
   double precision, allocatable :: x(:,:,:), x0(:,:,:), xfinal(:), xtemp(:)

   ! For reading parameters
   integer iconf(7)
   double precision conf(2)

   ! Space and time steps
   double precision dt, dt1, dt2, hx, hy, hz

   ! Current error and limit convergence
   double precision resLoc, result, epsilon

   ! Index variables
   integer i, j, k, l, p

   ! Time and step variables
   double precision t
   integer step

   ! Max step
   integer maxStep

   ! Variables for clock
   double precision time_init, time_final, elapsed_time

   ! Number of initial borders layers to
   ! avoid artifacts problem on corners
   integer nb_layers

   ! MPI variables
   integer, dimension(3) :: sizes, subsizes1, subsizes2, subsizes3, starts
   integer nproc, ndims, infompi, comm, comm3d
   integer, dimension(3) :: dims
   logical, dimension(3) :: periods
   logical, parameter :: reorganisation = .false.
   integer matrix_type_oxz, matrix_type_oxy, matrix_type_oyz
   integer, parameter :: nbvi=4
   integer, parameter :: S=1, E=2, N=3, W=4, Zd=5, Zu=6
   integer, dimension(6) :: neighBor
   integer xcell, ycell, zcell, size_tot_x, size_tot_y, size_tot_z
   integer, allocatable :: xs(:), ys(:), zs(:), xe(:), ye(:), ze(:)

   ! Physical parameters
   double precision temp1_init, temp2_init, k0

   ! Initial number of borders layers
   nb_layers = 1

   ! temp1_init: temperature init on borders
   temp1_init = 10.0

   ! temp2_init: temperature init inside
   temp2_init = -10.0

   ! Diffusivity coefficient
   k0 = 1

   ! MPI initialization
   call MPI_Init(infompi)
   comm = MPI_COMM_WORLD
   call MPI_Comm_size(comm, nproc, infompi)
   call MPI_Comm_rank(comm, me, infompi)

   ! Get input parameters
   if(me.eq.0) then
    call readParam(iconf, conf)
   endif

   ! Broadcast input parameters
   call MPI_Bcast(iconf, 7, MPI_INTEGER, 0, comm, infompi)
   call MPI_Bcast(conf, 2, MPI_DOUBLE_PRECISION, 0, comm, infompi)

   ! Assign input parameters to variables
   size_x    = iconf(1)
   size_y    = iconf(2)
   size_z    = iconf(3)
   x_domains = iconf(4)
   y_domains = iconf(5)
   z_domains = iconf(6)
   maxStep   = iconf(7)
   dt1       = conf(1)
   epsilon   = conf(2)

   ! Warning message if dimensions and number of processes don't match
   if((me.eq.0).and.(nproc.ne.(x_domains*y_domains*z_domains))) then
    write(*,*) 'Number of processes not equal to Number of subdomains'
   endif

   ! Various other variables
   size_x_glo = size_x+2
   size_y_glo = size_y+2
   size_z_glo = size_z+2
   hx = 1.0d0/dble(size_x_glo)
   hy = 1.0d0/dble(size_y_glo)
   hz = 1.0d0/dble(size_z_glo)
   dt2 = 0.125*(min(min(hx,hy),hz)**2)/k0
   size_tot_x = size_x+2*x_domains+2
   size_tot_y = size_y+2*y_domains+2
   size_tot_z = size_z+2*z_domains+2

   ! Take a right time step for convergence
   if(dt1.ge.dt2) then
    if(me.eq.0) then
     write(*,*)
     write(*,*) 'Time step too large in param file - Taking convergence criterion'
    endif
    dt = dt2
   else
    dt = dt1
   endif

   ! Allocate 3D contiguous arrays
   allocate(xfinal(1:size_x*size_y*size_z))
   allocate(x(0:size_tot_x-1,0:size_tot_y-1,0:size_tot_z-1))
   allocate(x0(0:size_tot_x-1,0:size_tot_y-1,0:size_tot_z-1))

   ! Allocate coordinates of processes
   allocate(xs(0:nproc-1))
   allocate(xe(0:nproc-1))
   allocate(ys(0:nproc-1))
   allocate(ye(0:nproc-1))
   allocate(zs(0:nproc-1))
   allocate(ze(0:nproc-1))

   ! Create 3D cartesian grid
   periods(:) = .false.

   ndims = 3
   dims(1) = x_domains
   dims(2) = y_domains
   dims(3) = z_domains

   call MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, &
                        reorganisation, comm3d, infompi)

   ! Identify neighbors
   NeighBor(:) = MPI_PROC_NULL

   ! Left/West and right/Est neighbors
   call MPI_Cart_shift(comm3d, 0, 1, NeighBor(W), NeighBor(E), infompi)

   ! Bottom/South and Upper/North neighbors
   call MPI_Cart_shift(comm3d, 1, 1, NeighBor(S), NeighBor(N), infompi)

   ! Zdown/South and Zup/North neighbors
   call MPI_Cart_shift(comm3d, 2, 1, NeighBor(Zd), NeighBor(Zu), infompi)

   ! Size of each cell
   xcell = (size_x/x_domains)
   ycell = (size_y/y_domains)
   zcell = (size_z/z_domains)

   ! Allocate subdomain
   allocate(xtemp(1:xcell*ycell*zcell))

   ! Compute xs, xe, ys, ye, zs, ze for each cell on the grid
   call processToMap(xs, xe, ys, ye, zs, ze, xcell, ycell, zcell, x_domains, y_domains, z_domains, nproc)

   ! Create matrix data types to communicate
   sizes(1) = size_tot_x
   sizes(2) = size_tot_y
   sizes(3) = size_tot_z

   starts(1) = 0
   starts(2) = 0
   starts(3) = 0

   ! Create matrix data type to communicate on vertical Oxz plane
   subsizes1(1) = xcell
   subsizes1(2) = 1
   subsizes1(3) = zcell

   call MPI_Type_create_subarray(3, sizes, subsizes1, starts, &
                                 MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
                                 matrix_type_oxz, infompi)
   call MPI_Type_commit(matrix_type_oxz, infompi)

   ! Create matrix data type to communicate on vertical Oyz plane
   subsizes2(1) = 1
   subsizes2(2) = ycell
   subsizes2(3) = zcell

   call MPI_Type_create_subarray(3, sizes, subsizes2, starts, &
                                 MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
                                 matrix_type_oyz, infompi)
   call MPI_Type_commit(matrix_type_oyz, infompi)

   ! Create matrix data type to communicate on vertical Oxy plane
   subsizes3(1) = xcell
   subsizes3(2) = ycell
   subsizes3(3) = 1

   call MPI_Type_create_subarray(3, sizes, subsizes3, starts, &
                                 MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
                                 matrix_type_oxy, infompi)
   call MPI_Type_commit(matrix_type_oxy, infompi)

   ! Initialize values
   call initValues(nb_layers, x0, size_tot_x, size_tot_y, size_tot_z, temp1_init, temp2_init)

   ! Update the boundaries
   call updateBound(x0, size_tot_x, size_tot_y, size_tot_z, NeighBor, comm3d, &
                    matrix_type_oxz, matrix_type_oxy, matrix_type_oyz, me, &
                    xs, ys, zs, xe, ye, ze, nproc)

   ! Initialize step and time
   step = 0
   t = 0.0

   ! Starting time
   time_init = MPI_Wtime()

   ! Main loop
10 continue

   ! Increment step and time
   step = step + 1
   t = t + dt

   ! Perform one step of the explicit scheme
   call computeNext(x0, x, size_tot_x, size_tot_y, size_tot_z, dt, hx, hy, hz, resLoc, &
                    me, xs, ys, zs, xe, ye, ze, nproc, k0)

   ! Update the partial solution along the interface
   call updateBound(x0, size_tot_x, size_tot_y, size_tot_z, NeighBor, comm3d, &
                    matrix_type_oxz, matrix_type_oxy, matrix_type_oyz, &
                    me, xs, ys, zs, xe, ye, ze, nproc)

   ! Sum reduction to get error
   call MPI_Allreduce(resLoc, result, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, infompi)

   ! Current error
   result = sqrt(result)

   ! Break conditions of main loop
   if(.not.((result.lt.epsilon).or.(step.gt.maxStep))) goto 10

   ! Gather all subdomains
   i = 1
   do 301, k=zs(me),ze(me)
    l = 1
    do 300, j=ys(me),ye(me)

     xtemp((l-1)*xcell+(i-1)*xcell*ycell+1:l*xcell+(i-1)*xcell*ycell) = x0(xs(me):xe(me),j,k)
     l = l+1

300 continue
    i = i+1
301 continue

   ! Perform gathering
   call MPI_Gather(xtemp, xcell*ycell*zcell, MPI_DOUBLE_PRECISION, xfinal, xcell*ycell*zcell, &
                   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, infompi)

   ! Ending time
   time_final = MPI_Wtime()
   ! Elapsed time
   elapsed_time = time_final - time_init

   ! Print results
   if(me.eq.0) then
    write(*,*)
    write(*,*) '  Time step = ',dt
    write(*,*)
    write(*,1001) epsilon,step
    write(*,*)
    write(*,*) '  Problem size = ',size_x*size_y*size_z
    write(*,*)
    write(*,1002) elapsed_time
    write(*,*)
    write(*,*) '  Computed solution in outputPar.dat file '
    write(*,*)

    ! Store solution into output file
    open(5,file='outputPar.dat',action='write',status='replace')
    do 745, j=1,size_y+2
     write(5,1000,advance='no') (temp1_init,i=1,size_x+1)
     write(5,999,advance='no') temp1_init
     write(5,*)
745 continue
    write(5,*)

    do 612, p=1,size_z
     write(5,1000) (temp1_init,i=1,size_x+2)
     do 400, i=1,y_domains
      do 500, j=1,ycell
       write(5,1000,advance='no') temp1_init
       do 600, k=0,x_domains-1

        write(5,1000,advance='no') (xfinal((j-1)*xcell+l+(y_domains-i)*z_domains*xcell*ycell*zcell+&
                                    k*y_domains*z_domains*xcell*ycell*zcell+(p-1)*xcell*ycell),l=1,xcell)
600    continue
       write(5,999,advance='no') temp1_init
       write(5,*)
500   continue
400  continue
     write(5,1000) (temp1_init,i=1,size_x+2)
     write(5,*)
612 continue

    do 746, j=1,size_y+2
     write(5,1000,advance='no') (temp1_init,i=1,size_x+1)
     write(5,999,advance='no') temp1_init
     write(5,*)
746 continue
    close(5)
   endif

   ! Free arrays
   deallocate(x)
   deallocate(x0)
   deallocate(xtemp)
   deallocate(xfinal)
   deallocate(xs)
   deallocate(xe)
   deallocate(ys)
   deallocate(ye)
   deallocate(zs)
   deallocate(ze)

   ! Free matrices type
   call MPI_Type_free(matrix_type_oxz,infompi)
   call MPI_Type_free(matrix_type_oyz,infompi)
   call MPI_Type_free(matrix_type_oxy,infompi)

   call MPI_Finalize(infompi)

   ! Formats available to display the computed values on the grid
999 format(1000(f15.11))
1000 format(1000(f15.11,1x))
1001 format('   Convergence = ',f11.9,' after ',i9,' steps ')
1002 format('   Wall Clock = ',f15.6)

   stop
end

!!explicit_Par ends
