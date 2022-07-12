program verify

use params

double precision, dimension(nx,ny) :: Ts,Tp
double precision :: time_s,time_p
integer :: it_s,it_p

integer :: i,j


open(343,file='../serial/info.dat',status='old',form='formatted')
read(343,*)
read(343,'(f20.3,i20)') time_s,it_s
close(343,status='keep')

open(343,file='../OpenMP/info.dat',status='old',form='formatted')
read(343,*)
read(343,'(f20.3,i20)') time_p,it_p
close(343,status='keep')


open(343,file='../serial/T.dat',status='old',form='formatted')
do j=ny,1,-1
  do i=1,nx
    read(343,'(es20.12)', advance='no') Ts(i,j)
  enddo
  read(343,*)
enddo
close(343,status='keep')

open(343,file='../OpenMP/T.dat',status='old',form='formatted')
do j=ny,1,-1
  do i=1,nx
    read(343,'(es20.12)', advance='no') Tp(i,j)
  enddo
  read(343,*)
enddo
close(343,status='keep')



write(*,'(3(a15))') '','Serial','Parallel'
write(*,'(a15,2(i15))') 'Iterations',it_s,it_p
write(*,'(a15,2(f15.3),a,f10.3)') 'Time [s]',time_s,time_p,'       speed-up  x',time_s/time_p
if(maxval(maxval(abs(Ts-Tp),2),1).lt.1e-15)then
 write(*,'(a)') 'Code verified, same output'
else
 write(*,'(a)') 'Code NOT verified, different output, maxerror: ',maxval(maxval(abs(Ts-Tp),2),1)
endif

end program verify
