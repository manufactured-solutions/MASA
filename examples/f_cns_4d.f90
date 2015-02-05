!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! MASA - Manufactured Analytical Solutions Abstraction Library
!!
!! Copyright (C) 2010,2011,2012,2013 The PECOS Development Team
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
!! f_cns_4d.f90: Fortran Compressible Navier Stokes 4D Example
!!
!! -------------------------------------------------------------------------
!! -------------------------------------------------------------------------

program main
  use masa
  use iso_c_binding

  implicit none

  integer :: err=0

  ! point to evaluate (just one for now)
  real(8) x, y, z, t

  ! source terms
  real(8) q_rho, q_ru, q_rv, q_rw, q_re

  ! solution
  real(8) rho, u, v, w, p

  ! init masa
  call masa_init('cns_4d_example','navierstokes_4d_compressible_powerlaw');
  call masa_init_param()
  call masa_sanity_check()

  ! some point
  x = 1.0d0
  y = 2.0d0
  z = 0.5d0
  t = 0.1d0

  ! evaluate source terms
  q_rho = masa_eval_4d_source_rho  (x,y,z,t)
  call test(q_rho)

  q_ru  = masa_eval_4d_source_rho_u(x,y,z,t)
  call test(q_ru)

  q_rv  = masa_eval_4d_source_rho_v(x,y,z,t)
  call test(q_rv)

  q_rw  = masa_eval_4d_source_rho_w(x,y,z,t)
  call test(q_rw)

  q_re  = masa_eval_4d_source_rho_e(x,y,z,t)
  call test(q_re)

  ! evaluate solution
  rho = masa_eval_4d_exact_rho(x,y,z,t)
  call test(rho)

  u   = masa_eval_4d_exact_u  (x,y,z,t)
  call test(u)

  v   = masa_eval_4d_exact_v  (x,y,z,t)
  call test(v)

  w   = masa_eval_4d_exact_w  (x,y,z,t)
  call test(w)

  p   = masa_eval_4d_exact_p  (x,y,z,t)
  call test(p)


  call exit(err)
end program main

subroutine test(in)
  implicit none

  real(8), intent(in) :: in
  real(8)             :: MASA_DEFAULT
  real(8)             :: MASA_UNINIT

  MASA_DEFAULT = -12345.67
  MASA_UNINIT  = -1.33d0

  if(in.eq.MASA_DEFAULT) then
     call exit(1)
  endif

  if(in.eq.MASA_UNINIT) then
     call exit(1)
  endif

end subroutine test
