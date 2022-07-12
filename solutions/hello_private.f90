program hello
integer i
!$ integer OMP_GET_THREAD_NUM

 i=-1

!$OMP PARALLEL private(i)
!$ i=OMP_GET_THREAD_NUM()
 print(*,*) 'Hello world', i
!$OMP END PARALLEL

end program hello
