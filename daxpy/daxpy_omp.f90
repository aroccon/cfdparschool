subroutine daxpy(n,a,x,y)
 double precision :: x(n),y(n),a
 integer :: n, i
!$omp target teams distribute
 do i=1,n
   y(i) = a*x(i)+y(i)
 enddo
end subroutine daxpy


program main
use omp_lib
implicit none
integer ::n=50000000,i
double precision :: a=3.0d0
double precision,allocatable :: x(:), y(:)
double precision :: t1,t2

allocate(x(n))
allocate(y(n))

call omp_set_default_device(1)
t1 = omp_get_wtime()
call daxpy(n,a,x,y)
t2 = omp_get_wtime()
write(*,*) 'Time elapsed:', t2-t1

end program
