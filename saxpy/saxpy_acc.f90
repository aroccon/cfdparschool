program main
integer ::n=1000,i
real:: a=3.0
real,allocatable :: x(:), y(:)

allocate(x(n))
allocate(y(n))

x=2.0d0
y=1.0d0

!!!!!!$acc kernels
!$acc parallel loop 
do i=1,n
    y(i)=a*x(i)+y(i)
enddo
!$acc end parallel
!!!!!!$acc end kernels

end program
