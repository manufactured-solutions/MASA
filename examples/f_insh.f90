!! -*-f90-*-
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
!! f_laplace.f90: Fortran Incompressible Navierstokes 3D Example
!!
!! $Id: 
!! -------------------------------------------------------------------------
!! -------------------------------------------------------------------------

program main
  use masa
  use iso_c_binding

  implicit none

  ! solutions
  real(8) :: u,v,w
  real(8) :: exact_u,exact_v,exact_w

  ! error handling
  integer ::  err=0 

  ! declarations
  real(8) :: x
  real(8) :: y
  real(8) :: z
  real(8) :: fsol

  ! problem size
  integer i,j,k
  integer ::  nx = 8
  integer ::  ny = 9
  integer ::  nz = 12  
  integer ::  lx=3
  integer ::  ly=1 
  integer ::  lz=2 
  real(8) ::  dx 
  real(8) ::  dy
  real(8) ::  dz

  ! initialize the problem
  dx = real(lx)/real(nx)
  dy = real(ly)/real(ny);
  dz = real(lz)/real(nz);

  ! initialize the problem
  call masa_init("incompressible homogeneous isotropic NS example","navierstokes_3d_incompressible")

  ! evaluate source terms (3D)
  do i=0, nx
     do j=0, ny
        do k=0, nz
         
           x = i*dx
           y = j*dy        
           z = k*dz
        
           ! evalulate source term
           u = masa_eval_3d_source_u   (x,y,z)
           v = masa_eval_3d_source_v   (x,y,z)
           w = masa_eval_3d_source_w   (x,y,z)
           
           !evaluate analytical term
           exact_u = masa_eval_3d_exact_u (x,y,z)
           exact_v = masa_eval_3d_exact_v (x,y,z)
           exact_w = masa_eval_3d_exact_w (x,y,z)
           
           call test(u)
           call test(v)
           call test(w)

           call test(exact_u)
           call test(exact_v)
           call test(exact_w)
           
        enddo
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
