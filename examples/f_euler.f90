!! -*-f90-*-
!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! MASA - Manufactured Analytical Solutions Abstraction Library
!!
!! Copyright (C) 2010 The PECOS Development Team
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
!! f_cns.f90: Fortran Euler Example
!!
!! $Id: 
!! -------------------------------------------------------------------------
!! -------------------------------------------------------------------------

program main
  use masa
  use iso_c_binding
  implicit none

  ! declarations
  real(8) :: value
  real(8) :: out 

  ! initialize the problem
  call masa_init('mytest','euler_1d')

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
  out = masa_eval_1d_source_u(value)
  write(6,*) "value is: ", out 

  call exit(0)
end program main
