module vecaddmod
  implicit none
 contains
  subroutine vecaddgpu( r, a, b, n )
   real, dimension(:) :: r, a, b
   integer :: n
   integer :: i
!$acc kernels loop copyin(a(1:n),b(1:n)) copyout(r(1:n))
   do i = 1, n
    r(i) = a(i) + b(i)
   enddo
  end subroutine
end module

program test2
  use vecaddmod
  implicit none
  integer :: n, i, errs, argcount
  real, dimension(:), allocatable :: a, b, r, e
  character*10 :: arg1
  argcount = command_argument_count()
  n = 1000000  ! default value
  if( argcount >= 1 )then
   call get_command_argument( 1, arg1 )
   read( arg1, '(i)' ) n
   if( n <= 0 ) n = 100000
  endif
  allocate( a(n), b(n), r(n), e(n) )
  do i = 1, n
   a(i) = i
   b(i) = 1000*i
  enddo
  ! compute on the GPU
  call vecaddgpu( r, a, b, n )
  ! compute on the host to compare
  do i = 1, n
   e(i) = a(i) + b(i)
  enddo
  ! compare results
  errs = 0
  do i = 1, n
   if( r(i) /= e(i) )then
     errs = errs + 1
   endif
  enddo
  print *, errs, ' errors found'
  if( errs ) call exit(errs)
end program
