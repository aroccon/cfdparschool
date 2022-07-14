subroutine daxpy(n,a,x,y)
  double precision :: x(n),y(n),a
  integer :: n,i
  do i=1,n
    y(i) = a*x(i)+y(i)
  enddo
end subroutine daxpy


program main
implicit none
integer :: n=50000000,i
double precision :: a
double precision, allocatable :: x(:), y(:)
double precision :: t1,t2

allocate(x(n),y(n))

x=1.0d0
y=2.0d0
a=3.0d0

call cpu_time(t1)
call daxpy(n,a,x,y)
call cpu_time(t2)
write(*,*) 'Time elapsed:', t2-t1

deallocate(x,y)

end program
