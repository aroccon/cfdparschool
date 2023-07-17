program pi
  implicit none
  integer, parameter :: blocks_per_image = 2**16
  integer, parameter :: block_size = 2**10
  real, dimension(block_size) :: x, y
  integer :: in_circle[*]
  integer :: i, n_circle, n_total
  real :: step, xfrom

  n_total = blocks_per_image * block_size * num_images()
  step = 1./real(num_images())
  xfrom = (this_image() - 1) * step
  in_circle = 0
  do i=1, blocks_per_image
     call random_number(x)
     call random_number(y)
     in_circle = in_circle + count((xfrom + step * x)** 2 + y**2 < 1.)
  end do
  sync all
  if (this_image() == 1) then
     n_circle = in_circle
     do i=2, num_images()
        n_circle = n_circle + in_circle[i]
     end do
     print *,"pi is approximately", 4.0d0*real(n_circle)/real(n_total), "exact", 4d0*atan(1.)
  end if
end program pi