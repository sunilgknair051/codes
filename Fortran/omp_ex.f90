! This example was taken from:
!
! http://en.wikipedia.org/wiki/OpenMP#Fortran_90_free_form
!
! Feel free to modify and distribute it!

program hello90
use omp_lib
implicit none
 
integer:: thread_id, nthreads

    !$omp parallel private(thread_id)

    thread_id = omp_get_thread_num()
    write (*,*) 'Hello World from thread', thread_id

    !$omp barrier
    if ( thread_id == 0 ) then
        nthreads = omp_get_num_threads()
        write (*,*) 'There are', nthreads, 'threads'
    end if
    
    !$omp end parallel

end program hello90