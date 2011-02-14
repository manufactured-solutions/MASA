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
!! f_euler1d.f90: 
!!
!! $Id: 
!! -------------------------------------------------------------------------
!! -------------------------------------------------------------------------

program main
  use masa
  implicit none

  real(8) :: u_0
  real(8) :: MASA_DEFAULT = -12345.67d0;

  ! start up problem
  call masa_init('mytest','euler_1d')

  ! problem auto_initalized, check correctly so
  u_0 = masa_get_param("u_0")
  if(u_0.eq.MASA_DEFAULT) then
     write(6,*) "MASA ERROR:: Variables not being auto initalized!"
     call exit(1)
  endif  

  ! purge them with fire
  call masa_purge_default_param()

  ! values should be set to default: checking
  u_0 = masa_get_param("u_0")  
  if(u_0.ne.MASA_DEFAULT) then
     write(6,*) "MASA ERROR:: Variables not being purged properly!"
     write(6,*) "Received: ", u_0
     write(6,*) "Instead : ", MASA_DEFAULT
     call exit(1)
  endif 

  ! steady as she goes (exit with success)
  call exit(0)

end program main
