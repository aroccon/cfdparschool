program heated_plate

use params

integer :: i,j,iteration

double precision, allocatable, dimension(:,:) :: T, T_old
double precision :: Lx=1.0d0, Ly=1.0d0
double precision :: dx,dy,dt
double precision :: error,change,delta
double precision :: tcpu0,tcpu1,t0,t1

integer :: omp_get_thread_num,omp_get_num_threads
double precision :: omp_get_wtime


allocate(T(0:nx+1,0:ny+1))
allocate(T_old(0:nx+1,0:ny+1))

! grid spacing
dx=Lx/dble(nx-1)
dy=Ly/dble(ny-1)

! time step, 0.4*(CFL)
dt=0.4*min((0.5*dx**2),(0.5*dy**2))

! initialize matrix
T=0.0d0

! set boundary conditions
!  column 0 and row ny+1 temperature is equal to zero
T(:,0)=0.0d0
T(nx+1,:)=0.0d0
!  column nx+1 and row zero linear temperature profile from 0 to 1
do j=1,ny
  T(nx+1,j)=dble(ny-j)/dble(ny-1)
enddo

do i=1,nx
  T(i,0)=dble(i-1)/dble(nx-1)
enddo

! initialize T_old
T_old=T


error=1.0d0
iteration=0
change=0.0d0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! START TIMED REGION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call cpu_time(tcpu0)
t0=omp_get_wtime()

!$OMP PARALLEL SHARED(error,iteration,change) private(i,j,delta)

do while (error.gt.errtol.and.error.lt.10.0d0)
!$OMP DO schedule(static) reduction(max:change)
  do j=1,ny
    do i=1,nx
      delta=dt*( (T_old(i+1,j)-2.0d0*T_old(i,j)+T_old(i-1,j))/dx**2 + &
                  (T_old(i,j+1)-2.0d0*T_old(i,j)+T_old(i,j-1))/dy**2  )
      T(i,j)=T_old(i,j)+delta
      change=max(delta,change)
    enddo
  enddo
!$OMP END DO
! implicit barrier (implies FLUSH) at end of parallel do region (unless you specify nowait clause)


!$OMP SINGLE
  error=change
  change=0.0d0
! just one thread updates iteration
  iteration=iteration+1
  ! write(*,*) iteration, error
!$OMP END SINGLE

!$OMP DO schedule(static)
  ! update T_old
  do j=1,ny
    do i=1,nx
      T_old(i,j)=T(i,j)
    enddo
  enddo
!$OMP END DO
enddo
!$OMP END PARALLEL

call cpu_time(tcpu1)
t1=omp_get_wtime()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! END TIMED REGION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! leave this separate parallel region: is outside timed region. Use it to check
! how many threads are actually being used (check set/export OMP_NUM_THREADS
! works correctly)
!$OMP PARALLEL
! write(*,*) 'Thread ID',omp_get_thread_num()
! OMP SINGLE is executed by only one thread (not necessarily the master thread)
!$OMP SINGLE
write(*,*) omp_get_num_threads()
!$OMP END SINGLE
!$OMP END PARALLEL

write(*,'(a,f5.2,a)') 'CPU time: ',tcpu1-tcpu0, ' seconds'
write(*,'(a,f5.2,a)') 'Wall time: ',t1-t0, ' seconds'

open(454,file='info.dat',status='replace',form='formatted')
write(454,'(2(a20))') 'Time [s]','iterations'
write(454,'(f20.3,i20)') t1-t0,iteration
close(454,status='keep')

call write_output(T,nx,ny)

deallocate(T,T_old)

end program heated_plate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine write_output(T,nx,ny)

integer :: nx,ny
integer :: i,j
double precision, dimension(0:nx+1,0:ny+1) :: T
character(len=200) :: namefile

write(namefile,'(a)') './T.dat'

open(456,file=trim(namefile),status='replace',form='formatted')

do j=ny,1,-1
  do i=1,nx
    write(456,'(es20.12)', advance='no') T(i,j)
  enddo
  write(456,*)
enddo

close(456,status='keep')

return
end
