program pic
 implicit none
 integer :: i
 double precision :: dx,pi,x
 integer :: N=100

dx=1.0d0/dble(N)
pi=0.0d0
!$OMP PARALLEL private(x) shared(pi)
!$OMP DO REDUCTION(+:pi)
do i=1,N
  !compute midpoint
  x=dx*(dble(i)-0.5d0)
  !compute integral
  pi=pi+4.0d0/(1.0d0+x**2)*dx
enddo
!$OMP END DO
!$OMP END PARALLEL

write(*,*) 'Computed value', pi

end program pic
