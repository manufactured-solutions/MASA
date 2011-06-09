program main
  use masa

  implicit none

  real(8) :: value
  real(8) :: out 

  call masa_init('mytest','euler_1d')
  call masa_init_param()

  call masa_list_mms()
  call masa_sanity_check()

  call masa_display_param

  print*,'getting original L' 
  value = masa_get_param("L")

  print*,'Original L = ',value

  call masa_set_param("L",3.1415d0)
  call masa_display_param

  print*,'getting updated L' 
  value = masa_get_param("L")

  print*,'Updated L = ',value

  out = masa_eval_1d_source_u(value)
  write(6,*) "value is: ", out 

end program
