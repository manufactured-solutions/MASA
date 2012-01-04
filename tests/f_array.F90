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
!! $Id: masa.f90 13217 2010-09-11 04:30:16Z nick $
!!
!! -------------------------------------------------------------------------
!! -------------------------------------------------------------------------

program main
  use masa
  implicit none

  integer(4)           :: n
  real(8),dimension(100) :: arr

  ! test display_array
  call masa_init('3rd test','radiation_integrated_intensity')
  !call masa_display_array()

  !call masa_get_array("vec_mean",n,arr)

end program
