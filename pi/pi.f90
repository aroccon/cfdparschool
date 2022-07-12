program pi
 implicit none
 integer :: i
 double precision :: dx,pi
 integer :: N=20

dx=1.0d0/dble(N)
pi=0.0d0
do i=1,N
  !compute midpoint
  x=dx*(dble(i)-0.5d0)
  !compute integral
  pi=pi+4.0d0/(1.0d0+x**2)*dx
enddo

end program pi
