program main
  use masa
  use iso_c_binding

  implicit none

  ! declarations
  real(8) :: value
  real(8) :: out 

  ! problem size


  ! initialize the problem
  call masa_init('mytest','navierstokes_2d_compressible')

  ! initialize the default parameters
  call masa_init_param()

  ! display selected manufactured solution
  call masa_list_mms()

  ! intialize the various parameters required for Euler 2D
  ! call the sanity check routine 
  ! (tests that all variables have been initialized)
  call masa_sanity_check()

  ! display all parameters and selected values
  call masa_display_param

  print*,'getting original L' 
  value = masa_get_param("L")

  print*,'Original L = ',value

  call masa_set_param("L",3.1415d0)
  call masa_display_param

  print*,'getting updated L' 
  value = masa_get_param("L")

  print*,'Updated L = ',value

  ! evaluate at a particular point
  out = masa_eval_2d_u_source(value,value)
  write(6,*) "value is: ", out 

  stop
end program main
