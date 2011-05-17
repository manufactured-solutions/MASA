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
!! f_cns.f90: Fortran Compressible Navier Stokes Example
!!
!! $Id: 
!! -------------------------------------------------------------------------
!! -------------------------------------------------------------------------

program main
  use masa
  use iso_c_binding

  implicit none

  ! solutions
  real(8) :: ufield,ufield2,ufield3
  real(8) :: vfield,vfield2,vfield3
  real(8) :: efield,efield2,efield3
  real(8) :: rho,rho2,rho3

  real(8) :: exact_u,exact_u2,exact_u3
  real(8) :: exact_v,exact_v2,exact_v3
  real(8) :: exact_p,exact_p2,exact_p3
  real(8) :: exact_rho,exact_rho2,exact_rho3

  ! error handling
  integer ::  err=0 

  ! declarations
  real(8) :: value, val
  real(8) :: out,out1=1
  real(8) :: x
  real(8) :: y
  real(8) :: fsol
  ! problem size
  integer i,j
  integer ::  nx = 71
  integer ::  ny = 93  
  integer ::  lx=3
  integer ::  ly=1 
  real(8) ::  dx 
  real(8) ::  dy

  ! initialize the problem
  dx = real(lx)/real(nx)
  dy = real(ly)/real(ny);

  ! initialize the problem
  call masa_init('mytest','navierstokes_2d_compressible')

  ! call the sanity check routine 
  ! (tests that all variables have been initialized)
  call masa_sanity_check()

  ! evaluate source terms (2D)
  do i=0, nx
     do j=0, ny
         
        y = j*dy        
        x = i*dx
        
        ! evalulate source terms
        ufield = masa_eval_2d_source_rho_u (x,y)
	vfield = masa_eval_2d_source_rho_v (x,y)
        efield = masa_eval_2d_source_rho_e (x,y)
        rho    = masa_eval_2d_source_rho   (x,y)

	!evaluate analytical terms
        exact_u = masa_eval_2d_exact_u     (x,y)
        exact_v = masa_eval_2d_exact_v     (x,y)
        exact_p = masa_eval_2d_exact_p     (x,y)
        exact_rho = masa_eval_2d_exact_rho (x,y)

        call test(ufield)
        call test(vfield)
        call test(efield)
        call test(rho)
        
        call test(exact_u)
        call test(exact_v)
        call test(exact_p)
        call test(exact_rho)

     enddo
  enddo

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
