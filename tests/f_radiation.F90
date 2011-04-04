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
!! f_radiation.f90: Fortran radiation test
!!
!! $Id: 
!! -------------------------------------------------------------------------
!! -------------------------------------------------------------------------

program main
  use masa

  implicit none

  ! threshold calculation
  real(8) :: thresh = 1e-15

  call masa_init("radiation","radiation_integrated_intensity");

  ! call the sanity check routine 
  ! (tests that all variables have been initialized)
  call masa_sanity_check()

  ! steady as she goes (exit with success)
  call exit(0)

end program main

