module mymod
 attributes(global) subroutine saxpy(n,a,x,y)
  real :: x(n),y(n),a
  integer :: n, i
  attribute(value) :: a, n
  i = threadIdx
  do i=1,n
    y(i) = a*x(i)+y(i)
  enddo
 end subroutine saxpy
end module mymod


program main
integer ::n=2**20,i
real:: a=3.0
real,allocatable :: x(:), y(:)

allocate(x(n))
allocate(y(n))

call cpu_time(t1)

call saxpy(n,a,x,y)

call cpu_time(t2)
write(*,*) 'Time elapsed:', t2-t1

end program
