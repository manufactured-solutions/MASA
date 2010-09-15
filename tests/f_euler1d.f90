program main
  use masa

  implicit none

  ! declarations
  real(8) :: value
  real(8) :: out,out1=1
  real(8) :: x = 1
  real(8) :: fsol
  ! problem size


  ! initialize the problem
  call masa_init('mytest','euler_1d')

  ! initialize the default parameters
  call masa_init_param()
  
  ! intialize the various parameters required for Euler 2D
  ! call the sanity check routine 
  ! (tests that all variables have been initialized)
  call masa_sanity_check()

  value = masa_get_param("L")
  call masa_set_param("L",3.1415d0)
  value = masa_get_param("L")
  write(6,*) value

! TODO: below here is funky

!!!  out1 = fsol(%VAL(x))
!!!  write(*,*) "out1 is ", out1
  ! evaluate at a particular point
!!!  out = masa_eval_1d_u_source(value)

  stop
end program main
