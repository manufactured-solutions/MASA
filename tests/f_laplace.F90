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
!! f_laplace.f90: Fortran laplace 2d test
!!
!! $Id: 
!! -------------------------------------------------------------------------
!! -------------------------------------------------------------------------

program main
  use masa
  implicit none

  ! threshold calculation
  real(8) :: thresh = 1d-14

  ! solutions
  real(8) :: field,field2,field3
  real(8) ::exact_phi,exact_phi2,exact_phi3

  ! variables 
  real(8) :: Lx
  real(8) :: Ly

  ! declarations

  real(8) :: value, val
  real(8) :: out,out1=1
  real(8) :: x
  real(8) :: y
  real(8) :: fsol

  ! problem size

  integer i,j,k
  integer ::  nx = 10
  integer ::  ny = 8  
  integer ::  l_x=2
  integer ::  l_y=1 
  real(8) ::  dx 
  real(8) ::  dy

  ! initialize the problem

  dx = 1.0d0*l_x/nx
  dy = 1.0d0*l_y/ny

  call masa_init('laplace regression test','laplace_2d')

  ! initialize the default parameters

  call masa_init_param()

  ! intialize the various parameters required for laplace 2D
  Lx = masa_get_param("Lx")
  Ly = masa_get_param("Ly")

  ! call the sanity check routine 
  ! (tests that all variables have been initialized)

  call masa_sanity_check()

  ! evaluate source terms (2D)
  do i=0, nx
     do j=0, ny

	x = i * dx
        y = j * dy
 
        ! evalulate source terms
        field     = masa_eval_2d_source_f (x,y)
        
        !evaluate analytical terms   
        exact_phi = masa_eval_2d_exact_phi(x,y)
        
        ! check against exact 
        field2     = 2d0*(Lx-x)**2d0 - 8d0*(Lx-x)*(Lx+x) + 2d0*(Lx+x)**2d0 
        field2     = field2 + 2d0*(Ly-y)**2d0 - 8d0*(Ly-y)*(Ly+y) + 2d0*(Ly+y)**2d0
        exact_phi2 = (Ly-y)**2d0 * (Ly+y)**2d0 + (Lx-x)**2d0 * (Lx+x)**2d0
 
#ifdef MASA_STRICT_REGRESSION

        field3     = abs(field-field2)        
        exact_phi3 = abs(exact_phi-exact_phi2)

#else

	field3     = abs(field-field2)/abs(field2)        
        exact_phi3 = abs(exact_phi-exact_phi2)/abs(exact_phi2)

#endif 

        ! source terms
        if(field3 .gt. thresh) then
           write(6,*) "FortMASA FATAL ERROR: laplace-2d"
           write(6,*) "f Field"
           write(6,*) "exceeded by: ", field3
           write(6,*) "masa:        ", field
           write(6,*) "maple:       ", field2
           write(6,*) "@ x,y:       ",x,y
           call exit(1)
        endif

        ! analytical terms now
        if(exact_phi3 .gt. thresh) then
           write(6,*) "FortMASA FATAL ERROR: laplace-2d"
           write(6,*) "phi an"
           call exit(1)
        endif

     enddo
  enddo

  ! steady as she goes (exit with success)

  call exit(0)

end program main
