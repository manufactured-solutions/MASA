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
!! f_euler2d.f90: Fortran Euler 2d test
!!
!! $Id: masa.f90 13217 2010-09-11 04:30:16Z nick $
!! -------------------------------------------------------------------------
!! -------------------------------------------------------------------------

program main
  use masa
  implicit none

  ! threshold calculation
  real(8) :: thresh = 1e-15

  ! solutions
  real(8) :: tfield,tfield2,tfield3
  real(8) :: t_an,t_an2,t_an3
  ! variables 
  real(8) :: A_x
  real(8) :: B_y
  real(8) :: C_z
  real(8) :: k_0

  ! declarations
  real(8) :: x
  real(8) :: y
  real(8) :: z
  real(8) :: fsol
  ! problem size
  integer i,j,k
  integer ::  nx = 11
  integer ::  ny = 21  
  integer ::  nz = 21  
  integer ::  lx=3
  integer ::  ly=1 
  integer ::  lz=2 
  real(8) ::  dx 
  real(8) ::  dy
  real(8) ::  dz

  ! external functions
  real(8) :: eval_1d_t_source
  real(8) :: eval_1d_t_an
  real(8) :: eval_2d_t_source
  real(8) :: eval_2d_t_an
  real(8) :: eval_3d_t_source
  real(8) :: eval_3d_t_an

  ! initialize the problem
  dx = real(lx)/real(nx)
  dy = real(ly)/real(ny);
  dz = real(lz)/real(nz);

  ! ---------------------------------------------------------
  ! evaluate source terms (1D)
  ! ---------------------------------------------------------

  call masa_init("temp-test-1d","heateq_1d_steady_const");

  ! initialize the default parameters
  call masa_init_param()

  ! intialize the various parameters required for Euler 2D
  A_x = masa_get_param("A_x");
  k_0 = masa_get_param("k_0");

  ! call the sanity check routine 
  ! (tests that all variables have been initialized)
  call masa_sanity_check()

  do i=0, nx
	
     x = i * dx

     tfield = masa_eval_1d_t_source  (x)
     t_an   = masa_eval_1d_t_an      (x)

     tfield2 = eval_1d_t_source(%val(x),%val(A_x),%val(k_0)) 
     t_an2   = eval_1d_t_an    (%val(x),%val(A_x))

     tfield3 = abs(tfield-tfield2)
     t_an3 = abs(t_an-t_an2)

     ! just need error checker
     if(tfield3 .gt. thresh) then
        write(6,*) "FortMASA FATAL ERROR: heat-1d"
        write(6,*) "T Field"
        write(6,*) "exceeded by: ", tfield3
        write(6,*) "masa:        ", tfield
        write(6,*) "maple:       ", tfield2
        write(6,*) "@ x:         ",x
        call exit(1)
     endif

     ! analytical terms now
     if(t_an3 .gt. thresh) then
        write(6,*) "FortMASA FATAL ERROR: heat-1d"
        write(6,*) "T an"
        call exit(1)
     endif     

  enddo

  ! ---------------------------------------------------------
  ! evaluate source terms (2D)
  ! ---------------------------------------------------------
  call masa_init("temp-test-2d","heateq_2d_steady_const");
  call masa_init_param()
  A_x = masa_get_param("A_x");
  B_y = masa_get_param("B_y");
  k_0 = masa_get_param("k_0");
  call masa_sanity_check()

  do i=0, nx
     do j=0, ny
        
	x = i * dx
        y = j * dy

        tfield = masa_eval_2d_t_source  (x,y)
        t_an   = masa_eval_2d_t_an      (x,y)
        
        tfield2 = eval_2d_t_source(%val(x),%val(y),%val(A_x),%val(B_y),%val(k_0)) 
        t_an2   = eval_2d_t_an    (%val(x),%val(y),%val(A_x),%val(B_y))
        
        tfield3 = abs(tfield-tfield2)
        t_an3 = abs(t_an-t_an2)
 
        if(tfield3 .gt. thresh) then
           write(6,*) "FortMASA FATAL ERROR: heat-2d"
           write(6,*) "T Field"
           write(6,*) "exceeded by: ", tfield3
           write(6,*) "masa:        ", tfield
           write(6,*) "maple:       ", tfield2
           write(6,*) "@ x,y:       ",x,y
           call exit(1)
        endif

        ! analytical terms now
        if(t_an3 .gt. thresh) then
           write(6,*) "FortMASA FATAL ERROR: heat-2d"
           write(6,*) "T an"
           write(6,*) "@ x,y:       ",x,y
           call exit(1)
        endif        

     enddo
  enddo

  ! ---------------------------------------------------------
  ! evaluate source terms (3D)
  ! ---------------------------------------------------------
  call masa_init("temp-test-3d","heateq_3d_steady_const");
  call masa_init_param()
  A_x = masa_get_param("A_x");
  B_y = masa_get_param("B_y");
  C_z = masa_get_param("C_z");
  k_0 = masa_get_param("k_0");
  call masa_sanity_check()

  do i=0, nx
     do j=0, ny
        do k=0, nz
        
           x = i * dx
           y = j * dy
           z = k * dz

           tfield = masa_eval_3d_t_source  (x,y,z)
           t_an   = masa_eval_3d_t_an      (x,y,z)

           tfield2 = eval_3d_t_source(%val(x),%val(y),%val(z),%val(A_x),%val(B_y),%val(C_z),%val(k_0)) 
           t_an2   = eval_3d_t_an    (%val(x),%val(y),%val(z),%val(A_x),%val(B_y),%val(C_z))

           tfield3 = abs(tfield-tfield2)
           t_an3 = abs(t_an-t_an2)

           if(tfield3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: heat-3d"
              write(6,*) "T Field"
              write(6,*) "exceeded by: ", tfield3
              write(6,*) "masa:        ", tfield
              write(6,*) "maple:       ", tfield2
              write(6,*) "@ x,y,z:     ",x,y
              call exit(1)
           endif

           ! analytical terms now
           if(t_an3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: heat-3d"
              write(6,*) "T an"
              write(6,*) "@ x,y,z:     ",x,y,z
              call exit(1)
           endif

        enddo
     enddo
  enddo
  ! steady as she goes (exit with success)
  call exit(0)

end program main
