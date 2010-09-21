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
  real(8) :: value, val
  real(8) :: out,out1=1
  real(8) :: x
  real(8) :: y
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

  enddo ! done with 1d

!   ! evaluate source terms (2D)
!   do i=0, nx
!      do j=0, ny
        
! 	x = i * dx
!         y = j * dy
               
!         ! evalulate source terms
!         ufield = masa_eval_2d_u_source  (x,y)
! 	vfield = masa_eval_2d_v_source  (x,y)
!         efield = masa_eval_2d_e_source  (x,y)
!         rho    = masa_eval_2d_rho_source(x,y)
	
! 	!evaluate analytical terms
! 	u_an = masa_eval_2d_u_an        (x,y)
! 	v_an = masa_eval_2d_v_an        (x,y)
! 	p_an = masa_eval_2d_p_an        (x,y)
! 	rho_an = masa_eval_2d_rho_an    (x,y)
		  
!         ! check against maple
!         ufield2 = eval_2d_u_source  (%val(x),%val(y),%val(u_0),%val(u_x),%val(u_y),%val(v_0),%val(v_x),%val(v_y), &
!              %val(rho_0),%val(rho_x),%val(rho_y),%val(p_0),%val(p_x),%val(p_y),%val(a_px),%val(a_py), &
!              %val(a_rhox),%val(a_rhoy),%val(a_ux),%val(a_uy),%val(a_vx),%val(a_vy),%val(L))

! 	vfield2 = eval_2d_v_source  (%val(x),%val(y),%val(u_0),%val(u_x),%val(u_y),%val(v_0),%val(v_x),%val(v_y), &
!              %val(rho_0),%val(rho_x),%val(rho_y),%val(p_0),%val(p_x),%val(p_y),%val(a_px),%val(a_py), &
!              %val(a_rhox),%val(a_rhoy),%val(a_ux),%val(a_uy),%val(a_vx),%val(a_vy),%val(L))

!         rho2    = eval_2d_rho_source(%val(x),%val(y),%val(u_0),%val(u_x),%val(u_y),%val(v_0),%val(v_x),%val(v_y), &
!              %val(rho_0),%val(rho_x),%val(rho_y),%val(p_0),%val(p_x),%val(p_y),%val(a_px),%val(a_py), &
!              %val(a_rhox),%val(a_rhoy),%val(a_ux),%val(a_uy),%val(a_vx),%val(a_vy),%val(L))
 
!         efield2 = eval_2d_e_source  (%val(x),%val(y),%val(u_0),%val(u_x),%val(u_y),%val(v_0),%val(v_x),%val(v_y), &
!              %val(rho_0),%val(rho_x),%val(rho_y),%val(p_0),%val(p_x),%val(p_y),%val(a_px),%val(a_py), &
!              %val(a_rhox),%val(a_rhoy),%val(a_ux),%val(a_uy),%val(a_vx),%val(a_vy),%val(Gamma),%val(mu),%val(L))

!         u_an2   = eval_2d_u_an  (%val(x),%val(y),%val(u_0),%val(u_x),%val(u_y),%val(a_ux),%val(a_uy),%val(L))
!         v_an2   = eval_2d_v_an  (%val(x),%val(y),%val(v_0),%val(v_x),%val(v_y),%val(a_vx),%val(a_vy),%val(L))
!         rho_an2 = eval_2d_rho_an(%val(x),%val(y),%val(rho_0),%val(rho_x),%val(rho_y),%val(a_rhox),%val(a_rhoy),%val(L))
!         p_an2   = eval_2d_p_an  (%val(x),%val(y),%val(p_0),%val(p_x),%val(p_y),%val(a_px),%val(a_py),%val(L))

!         ! need to add strict / non-strict regressions
!         ufield3 = abs(ufield-ufield2)
! 	vfield3 = abs(vfield-vfield2)
! 	efield3 = abs(efield-efield2)
! 	rho3    = abs(rho-rho2)
	
! 	u_an3   = abs(u_an-u_an2)
! 	v_an3   = abs(v_an-v_an2)
! 	rho_an3 = abs(rho_an-rho_an2)
! 	p_an3   = abs(p_an-p_an2)
 
!         ! just need error checker
!         if(ufield3 .gt. thresh) then
!            write(6,*) "FortMASA FATAL ERROR: euler-2d"
!            write(6,*) "U Field"
!            write(6,*) "exceeded by: ", ufield3
!            write(6,*) "masa:        ", ufield
!            write(6,*) "maple:       ", ufield2
!            write(6,*) "@ x,y:       ",x,y
!            call exit(1)
!         endif
        
!         if(vfield3 .gt. thresh) then
!            write(6,*) "FortMASA FATAL ERROR: euler-2d"
!            write(6,*) "V Field"
!            call exit(1)
!         endif
        
!         if(efield3 .gt. thresh) then
!            write(6,*) "FortMASA FATAL ERROR: euler-2d"
!            write(6,*) "E Field"
!            call exit(1)
!         endif

!         if(rho3 .gt. thresh) then
!            write(6,*) "FortMASA FATAL ERROR: euler-2d"
!            write(6,*) "Rho Field"
!            call exit(1)
!         endif

!         ! analytical terms now
!         if(u_an3 .gt. thresh) then
!            write(6,*) "FortMASA FATAL ERROR: euler-2d"
!            write(6,*) "U an"
!            call exit(1)
!         endif

!         if(v_an3 .gt. thresh) then
!            write(6,*) "FortMASA FATAL ERROR: euler-2d"
!            write(6,*) "V an"
!            call exit(1)
!         endif

!         if(p_an3 .gt. thresh) then
!            write(6,*) "FortMASA FATAL ERROR: euler-2d"
!            write(6,*) "P an"
!            call exit(1)
!         endif

!         if(rho_an3 .gt. thresh) then
!            write(6,*) "FortMASA FATAL ERROR: euler-2d"
!            write(6,*) "Rho an"
!            call exit(1)
!         endif             

!      enddo
!  enddo
  
  ! steady as she goes (exit with success)
  call exit(0)

end program main
