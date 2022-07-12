
program jacobi

integer, parameter:: nx=1000,ny=1000
integer:: i,j,iter=1
double precision :: max_err=1.e-4
double precision :: max_iter=10000
double precision :: err=2.0d0
double precision :: t1,t2
double precision, dimension(0:nx+1,0:ny+1) :: A, A_new

!bc conditions
do i=1,nx
  A(i,0)=dble(i-1)/dble(nx-1)
enddo
do j=1,ny
  A(nx+1,j)=dble(ny-j)/dble(ny-1)
enddo

call cpu_time(t1)

!$acc data copy(A), create(A_new)
do while ( err .ge. max_err .and. iter .le. max_iter)
  !$acc parallel loop collapse(2)
  do i=1,nx
    do j=1,ny
       A_new(i,j)=0.25*(A(i+1,j)+A(i-1,j)+A(i,j+1)+A(i,j-1))
     enddo
   enddo
   !$acc end parallel
   err=0.0d0
   !$acc parallel loop collapse(2)  reduction(max:err)
   do i=1,nx
     do j=1,ny
       err = max(abs(A_new(i,j) - A(i,j)),err)
       A(i,j) = A_new(i,j)
     enddo
   enddo
   !$acc end parallel
!print*,'error:',err
iter = iter+1
enddo
!$acc end data

call cpu_time(t2)
print*,'cpu time:',t2-t1

end program
