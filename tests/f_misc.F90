!! -*-f90-*-
!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! MASA - Manufactured Analytical Solutions Abstraction Library
!!
!! Copyright (C) 2010,2011 The PECOS Development Team
!!
!! This library is free software; you can redistribute it and/or
!! modify it under the terms of the Version 2.1 GNU Lesser General
!! Public License as published by the Free Software Foundation.
!!
!! This library is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!! Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public
!! License along with this library; if not, write to the Free Software
!! Foundation, Inc. 51 Franklin Street, Fifth Floor,
!! Boston, MA  02110-1301  USA
!!
!!-----------------------------------------------------------------------el-
!!
!! $Id: masa.f90 13217 2010-09-11 04:30:16Z nick $
!!
!! -------------------------------------------------------------------------
!! -------------------------------------------------------------------------

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

  ! test display_array
  call masa_init('3rd test','radiation_integrated_intensity')
  call masa_display_array()

end program
