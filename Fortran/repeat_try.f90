program repeat

 implicit none

 integer i,j,x,y,nx,ny,ndir
 parameter(nx=2,ny=2,ndir=8)
 real f(1:nx,1:ny,0:ndir),fn(1:nx,1:ny,0:ndir)

 DO j= 1,2
   call rain
 ENDDO

end program repeat

subroutine rain(nx)
  integer i,j,x,y,nx,ny,ndir
  parameter(nx=2,ny=2,ndir=8)
  real f(1:nx,1:ny,0:ndir),fn(1:nx,1:ny,0:ndir)
  DO x= 1,nx
    DO y =1,ny
      DO i = 0,ndir
        f(x,y,i) = x*y*i
        print *, f(x,y,i)
      END DO
    END DO
  END DO
END subroutine rain
