program main
    implicit none
    write (*,*) "Doing task A", this_image()
    call execute_command_line('')
    sync all
    write (*,*) "Doing task B", this_image()
    call execute_command_line('')
  end program main



  
