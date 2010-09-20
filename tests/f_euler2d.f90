! -----------------------------------------------------------------------bl-
! --------------------------------------------------------------------------
! 
!  MASA - Manufactured Analytical Solutions Abstraction Library
! 
!  Copyright (C) 2010 The PECOS Development Team
! 
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the Version 2.1 GNU Lesser General
!  Public License as published by the Free Software Foundation.
! 
!  This library is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!  Lesser General Public License for more details.
! 
!  You should have received a copy of the GNU Lesser General Public
!  License along with this library; if not, write to the Free Software
!  Foundation, Inc. 51 Franklin Street, Fifth Floor,
!  Boston, MA  02110-1301  USA
! 
! -----------------------------------------------------------------------el-
!
! f_euler2d.f90: Fortran Euler 2d test
!
! $Id: masa.f90 13217 2010-09-11 04:30:16Z nick $
! --------------------------------------------------------------------------
! --------------------------------------------------------------------------

program main
  use masa

  implicit none

  ! declarations
  real(8) :: value, val
  real(8) :: out,out1=1
  real(8) :: x = 1
  real(8) :: fsol
  ! problem size

  real(8) :: eval_2d_u_source

  ! initialize the problem
  call masa_init('mytest','euler_1d')

  ! initialize the default parameters
  call masa_init_param()
  
  ! intialize the various parameters required for Euler 2D
  ! call the sanity check routine 
  ! (tests that all variables have been initialized)
  call masa_sanity_check()

  val = eval_2d_u_source(val,val)
  write(6,*) val

  stop
end program main
