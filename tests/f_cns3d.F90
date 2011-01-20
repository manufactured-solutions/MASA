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
!! f_cns3d.f90: Fortran compressible navier stokes 3d test
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

  real(8) ::u_an,u_an2,u_an3
  real(8) ::v_an,v_an2,v_an3
  real(8) ::w_an,w_an2,w_an3
  real(8) ::p_an,p_an2,p_an3
  real(8) ::rho_an,rho_an2,rho_an3

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
  real(8) :: R
  real(8) :: k2
   
  ! declarations
  real(8) :: value, val
  real(8) :: x
  real(8) :: y
  real(8) :: z
  real(8) :: fsol
  ! problem size
  integer i,j,k
  integer ::  nx = 21
  integer ::  ny = 17
  integer ::  nz = 34  
  integer ::  lx=3
  integer ::  ly=1 
  integer ::  lz=3 
  real(8) ::  dx 
  real(8) ::  dy
  real(8) ::  dz

  ! external functions
  real(8) :: eval_3d_u_source
  real(8) :: eval_3d_v_source
  real(8) :: eval_3d_w_source
  real(8) :: eval_3d_e_source
  real(8) :: eval_3d_rho_source

  real(8) :: eval_3d_u_an
  real(8) :: eval_3d_v_an
  real(8) :: eval_3d_w_an
  real(8) :: eval_3d_p_an
  real(8) :: eval_3d_rho_an

  ! initialize the problem
  dx = real(lx)/real(nx)
  dy = real(ly)/real(ny)
  dz = real(lz)/real(nz)

  call masa_init("navier-stokes-test","navierstokes_3d_compressible");

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

  R     = masa_get_param("R");
  k2    = masa_get_param("k");

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
           ufield = masa_eval_3d_u_source  (x,y,z)
           vfield = masa_eval_3d_v_source  (x,y,z)
           wfield = masa_eval_3d_w_source  (x,y,z)
           efield = masa_eval_3d_e_source  (x,y,z)
           rho    = masa_eval_3d_rho_source(x,y,z)

           !evaluate analytical terms
           u_an = masa_eval_3d_u_an        (x,y,z)
           v_an = masa_eval_3d_v_an        (x,y,z)
           w_an = masa_eval_3d_w_an        (x,y,z)
           p_an = masa_eval_3d_p_an        (x,y,z)
           rho_an = masa_eval_3d_rho_an    (x,y,z)

           ! check against maple
           ufield2 = eval_3d_u_source  (%val(x),%val(y),%val(z),%val(u_0),%val(u_x),%val(u_y),%val(u_z),%val(v_0),&
                %val(v_x),%val(v_y),%val(v_z),%val(w_0),%val(w_x),%val(w_y),%val(w_z), &
                %val(rho_0),%val(rho_x),%val(rho_y),%val(rho_z),%val(p_0),%val(p_x),%val(p_y),%val(p_z), &
                %val(a_px),%val(a_py),%val(a_pz),%val(a_rhox),%val(a_rhoy),%val(a_rhoz),&
                %val(a_ux),%val(a_uy),%val(a_uz),%val(a_vx),%val(a_vy),%val(a_vz),%val(a_wx),%val(a_wy),%val(a_wz),&
                %val(mu),%val(L),%val(R),%val(k2))

           vfield2 = eval_3d_v_source  (%val(x),%val(y),%val(z),%val(u_0),%val(u_x),%val(u_y),%val(u_z),%val(v_0),&
                %val(v_x),%val(v_y),%val(v_z),%val(w_0),%val(w_x),%val(w_y),%val(w_z), &
                %val(rho_0),%val(rho_x),%val(rho_y),%val(rho_z),%val(p_0),%val(p_x),%val(p_y),%val(p_z), &
                %val(a_px),%val(a_py),%val(a_pz),%val(a_rhox),%val(a_rhoy),%val(a_rhoz),&
                %val(a_ux),%val(a_uy),%val(a_uz),%val(a_vx),%val(a_vy),%val(a_vz),%val(a_wx),%val(a_wy),%val(a_wz),&
                %val(mu),%val(L),%val(R),%val(k2))

           wfield2 = eval_3d_w_source  (%val(x),%val(y),%val(z),%val(u_0),%val(u_x),%val(u_y),%val(u_z),%val(v_0),&
                %val(v_x),%val(v_y),%val(v_z),%val(w_0),%val(w_x),%val(w_y),%val(w_z), &
                %val(rho_0),%val(rho_x),%val(rho_y),%val(rho_z),%val(p_0),%val(p_x),%val(p_y),%val(p_z), &
                %val(a_px),%val(a_py),%val(a_pz),%val(a_rhox),%val(a_rhoy),%val(a_rhoz),&
                %val(a_ux),%val(a_uy),%val(a_uz),%val(a_vx),%val(a_vy),%val(a_vz),%val(a_wx),%val(a_wy),%val(a_wz),&
                %val(mu),%val(L),%val(R),%val(k2))

           rho2    = eval_3d_rho_source(%val(x),%val(y),%val(z),%val(u_0),%val(u_x),%val(u_y),%val(u_z),%val(v_0),&
                %val(v_x),%val(v_y),%val(v_z),%val(w_0),%val(w_x),%val(w_y),%val(w_z), &
                %val(rho_0),%val(rho_x),%val(rho_y),%val(rho_z),%val(p_0),%val(p_x),%val(p_y),%val(p_z), &
                %val(a_px),%val(a_py),%val(a_pz),%val(a_rhox),%val(a_rhoy),%val(a_rhoz),&
                %val(a_ux),%val(a_uy),%val(a_uz),%val(a_vx),%val(a_vy),%val(a_vz),%val(a_wx),%val(a_wy),&
                %val(a_wz),%val(mu),%val(L),%val(R),%val(k2))

           efield2 = eval_3d_e_source  (%val(x),%val(y),%val(z),%val(u_0),%val(u_x),%val(u_y),%val(u_z),%val(v_0),&
                %val(v_x),%val(v_y), &
                %val(v_z),%val(w_0),%val(w_x),%val(w_y),%val(w_z), &
                %val(rho_0),%val(rho_x),%val(rho_y),%val(rho_z),%val(p_0),%val(p_x),%val(p_y),%val(p_z), &
                %val(a_px),%val(a_py),%val(a_pz),%val(a_rhox),%val(a_rhoy),%val(a_rhoz),&
                %val(a_ux),%val(a_uy),%val(a_uz),%val(a_vx),%val(a_vy),%val(a_vz),%val(a_wx),%val(a_wy),%val(a_wz),&
                %val(mu),%val(Gamma),%val(L),%val(R),%val(k2))

           u_an2   = eval_3d_u_an  (%val(x),%val(y),%val(z),%val(u_0),%val(u_x),%val(u_y),%val(u_z),%val(a_ux),&
                %val(a_uy),%val(a_uz),%val(L))

           v_an2   = eval_3d_v_an  (%val(x),%val(y),%val(z),%val(v_0),%val(v_x),%val(v_y),%val(v_z),%val(a_vx),&
                %val(a_vy),%val(a_vz),%val(L))

           w_an2   = eval_3d_w_an  (%val(x),%val(y),%val(z),%val(w_0),%val(w_x),%val(w_y),%val(w_z),%val(a_wx),&
                %val(a_wy),%val(a_wz),%val(L))

           rho_an2 = eval_3d_rho_an(%val(x),%val(y),%val(z),%val(rho_0),%val(rho_x),%val(rho_y),%val(rho_z),&
                %val(a_rhox),%val(a_rhoy),%val(a_rhoz),%val(L))

           p_an2   = eval_3d_p_an  (%val(x),%val(y),%val(z),%val(p_0),%val(p_x),%val(p_y),%val(p_z),%val(a_px),&
                %val(a_py),%val(a_pz),%val(L))


#ifdef MASA_STRICT_REGRESSION
           ufield3 = abs(ufield-ufield2)
           vfield3 = abs(vfield-vfield2)
           wfield3 = abs(wfield-wfield2)
           efield3 = abs(efield-efield2)
           rho3    = abs(rho-rho2)

           u_an3   = abs(u_an-u_an2)
           v_an3   = abs(v_an-v_an2)
           w_an3   = abs(w_an-w_an2)
           rho_an3 = abs(rho_an-rho_an2)
           p_an3   = abs(p_an-p_an2)

#else

           ufield3 = abs(ufield-ufield2)/abs(ufield2)
           vfield3 = abs(vfield-vfield2)/abs(vfield2)
           wfield3 = abs(wfield-wfield2)/abs(wfield2)
           efield3 = abs(efield-efield2)/abs(efield2)
           rho3    = abs(rho-rho2)/abs(rho2)

           u_an3   = abs(u_an-u_an2)/abs(u_an2)
           v_an3   = abs(v_an-v_an2)/abs(v_an2)
           w_an3   = abs(w_an-w_an2)/abs(w_an2)
           rho_an3 = abs(rho_an-rho_an2)/abs(rho_an2)
           p_an3   = abs(p_an-p_an2)/abs(p_an2)
        
#endif 
           ! just need error checker
           if(ufield3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: compressible navier stokes-3d"
              write(6,*) "U Field"
              write(6,*) "exceeded by: ", ufield3
              write(6,*) "masa:        ", ufield
              write(6,*) "maple:       ", ufield2
              write(6,*) "@ x,y,z:     ",x,y,z
              call exit(1)
           endif

           if(vfield3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: compressible navier stokes-3d"
              write(6,*) "V Field"
              call exit(1)
           endif

           if(wfield3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: compressible navier stokes-3d"
              write(6,*) "W Field"
              write(6,*) "exceeded by: ", wfield3
              write(6,*) "masa:        ", wfield
              write(6,*) "maple:       ", wfield2
              write(6,*) "@ x,y,z:     ",x,y,z             
              call exit(1)
           endif

           if(efield3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: compressible navier stokes-3d"
              write(6,*) "E Field"
              call exit(1)
           endif

           if(rho3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: compressible navier stokes-3d"
              write(6,*) "Rho Field"
              call exit(1)
           endif

           ! analytical terms now
           if(u_an3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: compressible navier stokes-3d"
              write(6,*) "U an"
              call exit(1)
           endif

           if(v_an3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: compressible navier stokes-3d"
              write(6,*) "V an"
              call exit(1)
           endif

           if(w_an3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: compressible navier stokes-3d"
              write(6,*) "W an"
              call exit(1)
           endif

           if(p_an3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: compressible navier stokes-3d"
              write(6,*) "P an"
              call exit(1)
           endif

           if(rho_an3 .gt. thresh) then
              write(6,*) "FortMASA FATAL ERROR: compressible navier stokes-3d"
              write(6,*) "Rho an"
              call exit(1)
           endif

        enddo
     enddo
  enddo
  ! steady as she goes (exit with success)
  call exit(0)

end program main
