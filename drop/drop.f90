!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finite volume solver with forward in time, centered in space discretization  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program drop

use data

integer :: i,j,it
character(len=200) :: fnameu,fnamerho,fnamep   !,command,fnamev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialization                                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! time step
dt=0.00125

! grid spacing and axes
dx=Lx/dble(nx)
dy=Ly/dble(ny)

allocate(x(nx+2))
allocate(y(ny+2))

do i=1,nx+2
  x(i)=dx*(i-1.5)
enddo

do j=1,ny+2
  y(j)=dy*(j-1.5)
enddo

! boundary conditions (tangential velocity at boundaries
utop=0.0d0
ubottom=0.0d0
vleft=0.0d0
vright=0.0d0

! velocity, pressure and density field
allocate(u(nx+1,ny+2))
allocate(ut(nx+1,ny+2))
allocate(v(nx+2,ny+1))
allocate(vt(nx+2,ny+1))
allocate(p(nx+2,nx+2))
allocate(rho(nx+2,ny+2))
allocate(rt(nx+2,ny+2))

allocate(uint(nx,ny))
allocate(vint(nx,ny))

! temp arrays
allocate(tmp1(nx+2,ny+2))
allocate(tmp2(nx+2,ny+2))
allocate(ptmp(nx+2,nx+2))

! viscosity
m0=0.01

! gravity vector
gx=0.0d0
gy=-100.0d0

! initialize density field
radius=0.15 ! drop radius
xc=0.5 ! x-position of drop center
yc=0.7  ! y-position of drop center

rho1=1.0d0
rho2=2.0d0

do j=1,ny+2
  do i=1,nx+2
    if ( (x(i)-xc)**2+(y(j)-yc)**2 .le. radius**2)then
      rho(i,j)=rho2
    else
      rho(i,j)=rho1
    endif
  enddo
enddo




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! time loop                                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do it=0,nstep
  if(mod(it,20).eq.0) write(*,*) 'Iteration ',it,' out of ',nstep
  ! tangential velocity at boundaries (ghost nodes)
  u(:,1)=2.0d0*ubottom-u(:,2)
  u(:,ny+2)=2.0d0*utop-u(:,ny+1)
  v(1,:)=2.0d0*vleft-v(2,:)
  v(nx+2,:)=2.0d0*vright-v(nx+1,:)

  ! temporary velocity
  do j=2,ny+1
    do i=2,nx
      ut(i,j)=u(i,j)+dt*(-0.25d0*( ((u(i+1,j)+u(i,j))**2-(u(i,j)+u(i-1,j))**2)/dx+ &
         ((u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j))-(u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1)))/dy)+ &
         m0/(0.5d0*(rho(i+1,j)+rho(i,j)))*((u(i+1,j)-2.0d0*u(i,j)+u(i-1,j))/dx**2+(u(i,j+1)-2.0d0*u(i,j)+u(i,j-1))/dy**2 )+gx)
    enddo
  enddo

  do j=2,ny
    do i=2,nx+1
      vt(i,j)=v(i,j)+dt*(-0.25d0*( ((u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j))-(u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j)))/dx+ &
         ((v(i,j+1)+v(i,j))**2-(v(i,j)+v(i,j-1))**2)/dy)+ &
         m0/(0.5d0*(rho(i,j+1)+rho(i,j)))*((v(i+1,j)-2.0d0*v(i,j)+v(i-1,j))/dx**2+(v(i,j+1)-2.0d0*v(i,j)+v(i,j-1))/dy**2 )+gy)
    enddo
  enddo

  ! compute source terms
  ! artificial density field for the computation of pressure (artificial values in the ghost cells)
  rt=rho
  rt(:,1)=1000.0d0
  rt(:,ny+2)=1000.0d0
  rt(1,:)=1000.0d0
  rt(nx+2,:)=1000.0d0

  do j=2,ny+1
    do i=2,nx+1
      tmp1(i,j)=(0.5d0/dt)*( (ut(i,j)-ut(i-1,j))/dx+(vt(i,j)-vt(i,j-1))/dy)
      tmp2(i,j)=1.0d0/( (1.0d0/dx)*( 1.0d0/(dx*(rt(i+1,j)+rt(i,j)))+1.0d0/(dx*(rt(i-1,j)+rt(i,j))) )+ &
                (1.0d0/dy)*(1.0d0/(dy*(rt(i,j+1)+rt(i,j)))+1.0d0/(dy*(rt(i,j-1)+rt(i,j))) ) );
    enddo
  enddo

  ! pressure solver (iterative)
  do k=1,maxit
    ptmp=p
    do j=2,ny+1
      do i=1,nx+1
        p(i,j)=(1.0d0-beta)*p(i,j)+beta*tmp2(i,j)*((1.0d0/dx)*( p(i+1,j)/(dx*(rt(i+1,j)+rt(i,j)))+ &
          p(i-1,j)/(dx*(rt(i-1,j)+rt(i,j))) )+(1.0d0/dy)*( p(i,j+1)/(dy*(rt(i,j+1)+rt(i,j)))+ &
          p(i,j-1)/(dy*(rt(i,j-1)+rt(i,j))) ) - tmp1(i,j));
      enddo
    enddo
    if(maxval(maxval(abs(ptmp-p),2),1).lt.errtol) exit
  enddo

  ! correct u velocity
  do j=2,ny+1
    do i=2,nx
      u(i,j)=ut(i,j)-dt*(2.0d0/dx)*(p(i+1,j)-p(i,j))/(rho(i+1,j)+rho(i,j));
    enddo
  enddo

  do j=2,ny
    do i=2,nx+1
      v(i,j)=vt(i,j)-dt*(2.0d0/dy)*(p(i,j+1)-p(i,j))/(rho(i,j+1)+rho(i,j))
    enddo
  enddo

  ! advect density field (central difference+diffusion)
  rt=rho
  do j=2,ny+1
    do i=2,nx+1
      rho(i,j)=rt(i,j)-(0.5d0*dt/dx)*(u(i,j)*(rt(i+1,j)+rt(i,j))-u(i-1,j)*(rt(i-1,j)+rt(i,j)) )- &
           (0.5d0* dt/dy)*(v(i,j)*(rt(i,j+1)+rt(i,j))-v(i,j-1)*(rt(i,j-1)+rt(i,j))  )+ &
           (m0*dt/dx**2)*(rt(i+1,j)-2.0d0*rt(i,j)+rt(i-1,j))+(m0*dt/dy**2)*(rt(i,j+1)-2.0d0*rt(i,j)+rt(i,j-1));
    enddo
  enddo


  ! periodically print output
  if(mod(it,10).eq.0)then
    ! interpolate velocity at cell center
    do i=1,nx
      uint(i,:)=0.5d0*(u(i+1,2:ny+1)+u(i,2:ny+1))
    enddo
    do j=1,ny
      vint(:,j)=0.5d0*(v(2:nx+1,j+1)+v(2:nx+1,j))
    enddo

    write(fnameu,'(a,i5.5,a)') './output/u_',it,'.dat'
    open(51,file=trim(fnameu),status='replace',form='formatted')
    ! write(fnamev,'(a,i5.5,a)') './output/v_',it,'.dat'
    ! open(52,file=trim(fnamev),status='replace',form='formatted')
    write(fnamerho,'(a,i5.5,a)') './output/rho_',it,'.dat'
    open(53,file=trim(fnamerho),status='replace',form='formatted')
    write(fnamep,'(a,i5.5,a)') './output/p_',it,'.dat'
    open(54,file=trim(fnamep),status='replace',form='formatted')
    ! write output gnuplot files
    do j=2,ny+1
      do i=2,nx+1
        write(51,'(4(es16.6))') x(i),y(j),uint(i-1,j-1),vint(i-1,j-1)
        write(53,'(3(es16.6))') x(i),y(j),rho(i,j)
        write(54,'(3(es16.6))') x(i),y(j),p(i,j)
      enddo
      ! write(51,*)
      ! write(52,*)
      write(53,*)
      write(54,*)
    enddo
    close(51,status='keep')
    ! close(52,status='keep')
    close(53,status='keep')
    close(54,status='keep')

    call write_vtk(it)

    ! if(it.eq.0)then
    !  write(command,'(a,a,a,a,a)') 'gnuplot -p -e "filenamerho=''',trim(fnamerho),'''" -e "filenameu=''',trim(fnameu),'''" plot.plt'
    ! else
    !  write(command,'(a,a,a,a,a)') 'gnuplot -e "filenamerho=''',trim(fnamerho),'''" -e "filenameu=''',trim(fnameu),'''" plot.plt'
    ! endif
    ! write(command,'(a,a,a,a,a)') 'gnuplot -e "filenamerho=''',trim(fnamerho),'''" -e "filenameu=''',trim(fnameu),'''" plot.plt'

    ! call system(command)
    ! call system('sleep 1')
    ! call system('open f.png')

  endif


enddo




deallocate(x,y)
deallocate(u,ut,v,vt,p,rho)
deallocate(tmp1,tmp2,ptmp,rt)
deallocate(uint,vint)

end program drop



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Additional subroutines                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write paraview                                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_vtk(iteration)

use data

integer :: iteration
integer :: i,j,k,nfields
integer :: nz=1

character(len=200) :: filename
character(len=80) :: buffer
character(len=10) :: lf
character(len=8) :: str1,str2,str3
character(len=16) :: str4

! end of line character
lf=achar(10)

! number of fields to be written
nfields=3

write(filename,'(a,i8.8,a)') './vtk/pview_',iteration,'.vtk'

open(666,file=trim(filename),status='replace',form='unformatted',access='stream',convert='big_endian')

! start writing vtk file
! write everything in single precision, saves storage, not needed double precision for visualization purposes

! write header
buffer='# vtk DataFile Version 3.0'//lf
write(666) trim(buffer)
buffer='Advection'//lf
write(666) trim(buffer)
buffer='BINARY'//lf
write(666) trim(buffer)
buffer='DATASET RECTILINEAR_GRID'//lf
write(666) trim(buffer)


!write grid
write(str1(1:8),'(i8)') nx
write(str2(1:8),'(i8)') ny
write(str3(1:8),'(i8)') nz
buffer='DIMENSIONS '//str1//str2//str3//lf
write(666) trim(buffer)
buffer='X_COORDINATES '//str1//'  float'//lf
write(666) trim(buffer)
do i=2,nx+1
 write(666) real(x(i))
enddo
buffer='Y_COORDINATES '//str2//'  float'//lf ;
write(666) trim(buffer)
do j=2,ny+1
 write(666) real(y(j))
enddo
buffer='Z_COORDINATES '//str3//'  float'//lf ;
write(666) trim(buffer)
do k=1,nz
 write(666) 0.0
enddo

! write content (data format)
write(str4(1:16),'(i16)') nx*ny*nz
buffer='POINT_DATA '//str4//lf
write(666) trim(buffer)
write(str1(1:8),'(i8)') nfields
buffer='FIELD FieldData '//str1//lf
write(666) trim(buffer)


! write velocity field (vector, 2D)
write(str4(1:16),'(i16)') nx*ny*nz
buffer='U 3 '//str4//' float'//lf
write(666) trim(buffer)
do k=1,nz
 do j=1,ny
  do i=1,nx
   write(666) real(uint(i,j)),real(vint(i,j)),0.0
  enddo
 enddo
enddo

! write rho field (scalar)
write(str4(1:16),'(i16)') nx*ny
buffer='RHO 1 '//str4//' float'//lf
write(666) trim(buffer)
do k=1,nz
 do j=2,ny+1
  do i=2,nx+1
   write(666) real(rho(i,j))
  enddo
 enddo
enddo

! write p field (scalar)
write(str4(1:16),'(i16)') nx*ny
buffer='P 1 '//str4//' float'//lf
write(666) trim(buffer)
do k=1,nz
 do j=2,ny+1
  do i=2,nx+1
   write(666) real(p(i,j))
  enddo
 enddo
enddo

buffer=lf
write(666) trim(buffer)
close(666,status='keep')



return
end
