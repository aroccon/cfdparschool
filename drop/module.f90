module data
 double precision :: Lx=1.0d0, Ly=1.0d0
 double precision :: dx,dy,dt
 double precision, allocatable, dimension(:,:) :: u,v,p,rho,ut,vt,rt,tmp1,tmp2,ptmp,uint,vint
 double precision, allocatable, dimension(:) :: x,y
 double precision :: xc,yc,radius
 double precision :: rho1,rho2,m0
 double precision :: gx,gy
 double precision :: utop,ubottom,vleft,vright
 double precision :: errtol=1.0e-3, beta=1.2
 integer :: nx=128, ny=128
 integer :: nstep=200,maxit=200
end module data
