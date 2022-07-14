module mymod
contains
 attributes(global) subroutine daxpy(n,a,x,y)
  double precision :: x(n),y(n),a
  integer :: n,i
  attributes(value) :: a,n
  n = size(x)
  i = blockDim%x * (blockIdx%x - 1) + threadIdx%x
  if (i <= n) y(i) = y(i) + a*x(i)*x(i)
end subroutine daxpy
end module mymod


program main
use cudafor
use mymod
integer ::n=50000000,i,istat
double precision :: a
double precision, allocatable :: x(:),y(:)
double precision, device :: x_d(n), y_d(n)
double precision :: t1,t2

allocate(x(n),y(n))

x=1.0d0
y=2.0d0
a=3.d0

call cpu_time(t1)
x_d=x
y_d=y
call daxpy<<<256, 256>>>(n,a,x_d,y_d)
y=y_d
call cpu_time(t2)
write(*,*) 'Time elapsed (CUDA):', t2-t1

deallocate(x,y)

end program
