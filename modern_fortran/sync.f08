program main
    implicit none
    write (*,*) "Doing task A", this_image()
    sync all
    write (*,*) "Doing task B", this_image()
  end program main


