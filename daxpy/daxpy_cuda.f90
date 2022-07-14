module mymod
contains
 attributes(global) subroutine daxpy(n,a,x,y)
  double precision :: x(n),y(n),a
  integer :: n,i
  attributes(value) :: a,n
  n = size(x)
  i = blockDim%x * (blockIdx%x - 1) + threadIdx%x
  if (i <= n) y(i) = y(i) + a*x(i)*x(i)*x(i) + sqrt(x(i))
end subroutine daxpy
end module mymod


program main
use cudafor
use mymod
integer ::n=50000000,i,istat
double precision :: a
double precision, device :: x_d(n), y_d(n)
double precision :: t1,t2

x_d=1.0
y_d=2.0
a=3.0d0

call cpu_time(t1)
call daxpy<<<4096, 256>>>(2**25,a,x_d,y_d)
istat = cudaDeviceSynchronize()
call cpu_time(t2)
write(*,*) 'Time elapsed (CUDA):', t2-t1

end program
