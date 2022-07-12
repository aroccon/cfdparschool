program hello
integer i
!$ integer OMP_GET_THREAD_NUM

 i=-1

!$OMP PARALLEL
!$OMP CRITICAL
!$ i=OMP_GET_THREAD_NUM()
 print(*,*) 'Hello world', i
!$OMP END CRITICAL
!$OMP END PARALLEL

end program hello
