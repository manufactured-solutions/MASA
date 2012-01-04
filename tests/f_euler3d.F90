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
!! f_euler3d.f90: Fortran Euler 3d test
!!
!! $Id: masa.f90 13217 2010-09-11 04:30:16Z nick $
!! -------------------------------------------------------------------------
!! -------------------------------------------------------------------------

program main
  use euler_source_interface
  use masa
  implicit none

  ! threshold calculation
  real(8) :: thresh = 1e-15

  ! solutions
  real(8) ::ufield,ufield2,ufield3
  real(8) ::vfield,vfield2,vfield3
  real(8) ::wfield,wfield2,wfield3
  real(8) ::efield,efield2,efield3
  real(8) ::rho,rho2,rho3

  real(8) ::exact_u,exact_u2,exact_u3
  real(8) ::exact_v,exact_v2,exact_v3
  real(8) ::exact_w,exact_w2,exact_w3
  real(8) ::exact_p,exact_p2,exact_p3
  real(8) ::exact_rho,exact_rho2,exact_rho3

  ! variables 
  real(8) :: u_0
  real(8) :: u_x
  real(8) :: u_y
  real(8) :: u_z
  real(8) :: v_0
  real(8) :: v_x
  real(8) :: v_y
  real(8) :: v_z
  real(8) :: w_0
  real(8) :: w_x
  real(8) :: w_y
  real(8) :: w_z
  real(8) :: rho_0
  real(8) :: rho_x
  real(8) :: rho_y
  real(8) :: rho_z
  real(8) :: p_0
  real(8) :: p_x
  real(8) :: p_y
  real(8) :: p_z
  real(8) :: a_px
  real(8) :: a_py
  real(8) :: a_pz
  real(8) :: a_rhox
  real(8) :: a_rhoy
  real(8) :: a_rhoz
  real(8) :: a_ux
  real(8) :: a_uy
  real(8) :: a_uz
  real(8) :: a_vx
  real(8) :: a_vy
  real(8) :: a_vz
  real(8) :: a_wx
  real(8) :: a_wy
  real(8) :: a_wz
  real(8) :: mu
  real(8) :: Gamma
  real(8) :: L    

  ! declarations
  real(8) :: value, val
  real(8) :: out,out1=1
  real(8) :: x
  real(8) :: y
  real(8) :: z
  real(8) :: fsol
  ! problem size
  integer i,j,k
  integer ::  nx = 71
  integer ::  ny = 17
  integer ::  nz = 34  
  integer ::  lx=3
  integer ::  ly=1 
  integer ::  lz=3 
  real(8) ::  dx 
  real(8) ::  dy
  real(8) ::  dz

  ! initialize the problem
  dx = real(lx)/real(nx)
  dy = real(ly)/real(ny)
  dz = real(lz)/real(nz)

  call masa_init('mytest','euler_3d')

  ! initialize the default parameters
  call masa_init_param()

  ! intialize the various parameters required for Euler 3D
  u_0 = masa_get_param("u_0");
  u_x = masa_get_param("u_x");
  u_y = masa_get_param("u_y");
  u_z = masa_get_param("u_z");

  v_0 = masa_get_param("v_0");
  v_x = masa_get_param("v_x");
  v_y = masa_get_param("v_y");
  v_z = masa_get_param("v_z");

  w_0 = masa_get_param("w_0");
  w_x = masa_get_param("w_x");
  w_y = masa_get_param("w_y");
  w_z = masa_get_param("w_z");

  rho_0 = masa_get_param("rho_0");
  rho_x = masa_get_param("rho_x");
  rho_y = masa_get_param("rho_y");
  rho_z = masa_get_param("rho_z");

  p_0 = masa_get_param("p_0");
  p_x = masa_get_param("p_x");
  p_y = masa_get_param("p_y");
  p_z = masa_get_param("p_z");

  a_px = masa_get_param("a_px");
  a_py = masa_get_param("a_py");
  a_pz = masa_get_param("a_pz");

  a_rhox = masa_get_param("a_rhox");
  a_rhoy = masa_get_param("a_rhoy");
  a_rhoz = masa_get_param("a_rhoz");

  a_ux = masa_get_param("a_ux");
  a_uy = masa_get_param("a_uy");
  a_uz = masa_get_param("a_uz");

  a_vx = masa_get_param("a_vx");
  a_vy = masa_get_param("a_vy");
  a_vz = masa_get_param("a_vz");

  a_wx = masa_get_param("a_wx");
  a_wy = masa_get_param("a_wy");
  a_wz = masa_get_param("a_wz");

  Gamma = masa_get_param("Gamma");
  mu    = masa_get_param("mu");
  L     = masa_get_param("L");

  ! call the sanity check routine 
  call masa_sanity_check()

  ! evaluate source terms (3D)
  do i=0, nx
     do j=0, ny
        do k=0, nz
           x = i * dx
           y = j * dy
           z = k * dz    

           ! evalulate source terms
           ufield = masa_eval_3d_source_rho_u  (x,y,z)
           vfield = masa_eval_3d_source_rho_v  (x,y,z)
           wfield = masa_eval_3d_source_rho_w  (x,y,z)
           efield = masa_eval_3d_source_rho_e  (x,y,z)
           rho    = masa_eval_3d_source_rho(x,y,z)

           !evaluate analytical terms
           exact_u = masa_eval_3d_exact_u        (x,y,z)
           exact_v = masa_eval_3d_exact_v        (x,y,z)
           exact_w = masa_eval_3d_exact_w        (x,y,z)
           exact_p = masa_eval_3d_exact_p        (x,y,z)
           exact_rho = masa_eval_3d_exact_rho    (x,y,z)

           ! need to hack this out
           ! check against maple
           ufield2 = eval_3d_source_u(x,y,z,u_0,u_x,u_y,&
                u_z,v_0,v_x,v_y, &
                v_z,w_0,w_x,w_y,w_z, &
                rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z, &
                a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,&
                a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy, &
                a_wz,L)

           vfield2 = eval_3d_source_v(x,y,z,u_0,u_x,u_y,&
                u_z,v_0,v_x,v_y, &
                v_z,w_0,w_x,w_y,w_z, &
                rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z, &
                a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,&
                a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,&
                a_wz,L)

           wfield2 = eval_3d_source_w(x,y,z,u_0,u_x,u_y,u_z,&
                v_0,v_x,v_y, &
                v_z,w_0,w_x,w_y,w_z, &
                rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z, &
                a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,&
                a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy, &
                a_wz,L)

           rho2    = eval_3d_source_rho(x,y,z,u_0,u_x,u_y,u_z,&
                v_0,v_x,v_y, &
                v_z,w_0,w_x,w_y,w_z, &
                rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z, &
                a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,&
                a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,&
                a_wz,mu,L)

           efield2 = eval_3d_source_e(x,y,z,u_0,u_x,u_y,u_z,&
                v_0,v_x,v_y, &
                v_z,w_0,w_x,w_y,w_z, &
                rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z, &
                a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,&
                a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,&
                mu,Gamma,L)

           exact_u2   = eval_3d_exact_u(x,y,z,u_0,u_x,u_y,u_z,a_ux,&
                a_uy,a_uz,L)
           
           exact_v2   = eval_3d_exact_v(x,y,z,v_0,v_x,v_y,v_z,a_vx,&
                a_vy,a_vz,L)

           exact_w2   = eval_3d_exact_w(x,y,z,w_0,w_x,w_y,w_z,a_wx,&
                a_wy,a_wz,L)

           exact_rho2 = eval_3d_exact_rho(x,y,z,rho_0,rho_x,rho_y,rho_z,&
                a_rhox,a_rhoy,a_rhoz,L)

           exact_p2   = eval_3d_exact_p(x,y,z,p_0,p_x,p_y,p_z,a_px,&
                a_py,a_pz,L)
           
#ifdef MASA_STRICT_REGRESSION
           ufield3 = abs(ufield-ufield2)
           vfield3 = abs(vfield-vfield2)
           wfield3 = abs(wfield-wfield2)
           efield3 = abs(efield-efield2)
           rho3    = abs(rho-rho2)

           exact_u3   = abs(exact_u-exact_u2)
           exact_v3   = abs(exact_v-exact_v2)
           exact_w3   = abs(exact_w-exact_w2)
           exact_rho3 = abs(exact_rho-exact_rho2)
           exact_p3   = abs(exact_p-exact_p2)

#else

           ufield3 = abs(ufield-ufield2)/abs(ufield2)
           vfield3 = abs(vfield-vfield2)/abs(vfield2)
           wfield3 = abs(wfield-wfield2)/abs(wfield2)
           efield3 = abs(efield-efield2)/abs(efield2)
           rho3    = abs(rho-rho2)/abs(rho2)

           exact_u3   = abs(exact_u-exact_u2)/abs(exact_u2)
           exact_v3   = abs(exact_v-exact_v2)/abs(exact_v2)
           exact_w3   = abs(exact_w-exact_w2)/abs(exact_w2)
           exact_rho3 = abs(exact_rho-exact_rho2)/abs(exact_rho2)
           exact_p3   = abs(exact_p-exact_p2)/abs(exact_p2)
        
#endif 

           ! just need error checker
           if(ufield3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: euler-3d"
              write(6,*) "U Field"
              write(6,*) "exceeded by: ", ufield3
              write(6,*) "masa:        ", ufield
              write(6,*) "maple:       ", ufield2
              write(6,*) "@ x,y,z:     ",x,y,z
              call exit(1)
           endif

           if(vfield3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: euler-3d"
              write(6,*) "V Field"
              call exit(1)
           endif

           if(wfield3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: euler-3d"
              write(6,*) "W Field"
              write(6,*) "exceeded by: ", wfield3
              write(6,*) "masa:        ", wfield
              write(6,*) "maple:       ", wfield2
              write(6,*) "@ x,y,z:     ",x,y,z             
              call exit(1)
           endif

           if(efield3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: euler-3d"
              write(6,*) "E Field"
              write(6,*) "exceeded by: ", efield3
              write(6,*) "masa:        ", efield
              write(6,*) "maple:       ", efield2
              write(6,*) "@ x,y,z:     ",x,y,z
              call exit(1)
           endif

           if(rho3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: euler-3d"
              write(6,*) "Rho Field"
              call exit(1)
           endif

           ! analytical terms now
           if(exact_u3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: euler-3d"
              write(6,*) "U an"
              call exit(1)
           endif

           if(exact_v3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: euler-3d"
              write(6,*) "V an"
              call exit(1)
           endif

           if(exact_w3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: euler-3d"
              write(6,*) "W an"
              call exit(1)
           endif

           if(exact_p3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: euler-3d"
              write(6,*) "P an"
              call exit(1)
           endif

           if(exact_rho3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: euler-3d"
              write(6,*) "Rho an"
              call exit(1)
           endif

        enddo
     enddo
  enddo
  ! steady as she goes (exit with success)
  call exit(0)

end program main
