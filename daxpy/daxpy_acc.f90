subroutine daxpy(n,a,x,y)
  double precision :: x(n),y(n),a
  integer :: n, i
!$acc data copyin(x,y,a) copyout(y)
!$acc kernels
  do i=1,n
    y(i) = a*x(i)+y(i)
  enddo
!$acc end kernels
!$acc end data
end subroutine daxpy


program main
use openacc
implicit none
integer ::n=50000000,i
double precision :: a=3.0d0
double precision,allocatable :: x(:), y(:)
double precision :: t1,t2

allocate(x(n))
allocate(y(n))

call acc_set_device_num(1,acc_device_nvidia)

call cpu_time(t1)
call daxpy(n,a,x,y)
call cpu_time(t2)
write(*,*) 'Time elapsed:', t2-t1

end program
