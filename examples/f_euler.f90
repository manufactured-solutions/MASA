!! -*-f90-*-
!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! MASA - Manufactured Analytical Solutions Abstraction Library
!!
!! Copyright (C) 2010,2011,2012 The PECOS Development Team
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
!! f_cns.f90: Fortran 1D Euler Example
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

  ! call the sanity check routine 
  ! (tests that all variables have been initialized)
  call masa_sanity_check()

  ! get original parameter value for L
  value = masa_get_param("L")

  ! set parameter value for L to be 3.14150000
  call masa_set_param("L",3.1415d0)

  ! check value is indeed 3.1415000 in masa
  value = masa_get_param("L")
  if(value.ne.3.1415d0) then
     write(6,*) "Incorrect value for L!"
     call exit(1)
  endif
  
  ! evaluate the source terms at a particular point
  out = masa_eval_1d_source_rho_u(value)

  call exit(0)
end program main
