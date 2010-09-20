program main
  use masa
  implicit none

  real(8) :: value = 0

  ! initialize the problem
  call masa_init('mytest','euler_1d')
  ! initialize the default parameters
  call masa_init_param()

  call masa_init('2nd test','euler_3d')
  call masa_init_param()
  
  ! test get_param function
  value = masa_get_param("L")

  if(value .ne. 3.02d0) then
     write(6,*) "FortMASA REGRESSION FAILURE: get param failure"
     write(6,*) "Exiting"
     call exit(1)
  endif 

  ! test set_param function
  call masa_set_param("L",3.1415d0)
  value = masa_get_param("L")

  if(value .ne. 3.1415d0) then
     write(6,*) "FortMASA REGRESSION FAILURE: set param failure"
     write(6,*) "Exiting"
     call exit(1)
  endif 

  ! test select_mms
  call masa_select_mms("mytest")
  value = masa_get_param("L") ! this should no longer be 3.14

  if(value .eq. 3.1415) then
     write(6,*) "FortMASA REGRESSION FAILURE: masa_switch failure"
     write(6,*) "Exiting"
     call exit(1)
  endif 

end program
