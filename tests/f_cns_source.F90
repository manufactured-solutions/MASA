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
!! $Author: nick $
!! $Id: f_euler_source.F90 16809 2011-01-20 00:21:14Z nick $
!!--------------------------------------------------------------------------
!!--------------------------------------------------------------------------

module cns_source_interface

  interface
     real (c_double) function eval_3d_u_source(x,y,z,&
          u_0,u_x,u_y,u_z,&
          v_0,v_x,v_y,v_z,&
          w_0,w_x,w_y,w_z,&
          rho_0,rho_x,rho_y,rho_z,&
          p_0,p_x,p_y,p_z, &
          a_px,a_py,a_pz,&
          a_rhox,a_rhoy,a_rhoz,&
          a_ux,a_uy,a_uz,&
          a_vx,a_vy,a_vz,&
          a_wx,a_wy,a_wz,&
          mu,L,R,k) bind (C) 

       use iso_c_binding
       implicit none

       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z
       real (c_double), value :: u_0
       real (c_double), value :: u_x
       real (c_double), value :: u_y
       real (c_double), value :: u_z
       real (c_double), value :: v_0
       real (c_double), value :: v_x
       real (c_double), value :: v_y
       real (c_double), value :: v_z
       real (c_double), value :: w_0
       real (c_double), value :: w_x
       real (c_double), value :: w_y
       real (c_double), value :: w_z
       real (c_double), value :: rho_0
       real (c_double), value :: rho_x
       real (c_double), value :: rho_y
       real (c_double), value :: rho_z
       real (c_double), value :: p_0
       real (c_double), value :: p_x
       real (c_double), value :: p_y
       real (c_double), value :: p_z
       real (c_double), value :: a_px
       real (c_double), value :: a_py
       real (c_double), value :: a_pz
       real (c_double), value :: a_rhox
       real (c_double), value :: a_rhoy
       real (c_double), value :: a_rhoz
       real (c_double), value :: a_ux
       real (c_double), value :: a_uy
       real (c_double), value :: a_uz
       real (c_double), value :: a_vx
       real (c_double), value :: a_vy
       real (c_double), value :: a_vz
       real (c_double), value :: a_wx
       real (c_double), value :: a_wy
       real (c_double), value :: a_wz
       real (c_double), value :: mu
       real (c_double), value :: L
       real (c_double), value :: R
       real (c_double), value :: k
       
     end function eval_3d_u_source

     real (c_double) function eval_3d_v_source(x,y,z,&
          u_0,u_x,u_y,u_z,&
          v_0,v_x,v_y,v_z,&
          w_0,w_x,w_y,w_z,&
          rho_0,rho_x,rho_y,rho_z,&
          p_0,p_x,p_y,p_z, &
          a_px,a_py,a_pz,&
          a_rhox,a_rhoy,a_rhoz,&
          a_ux,a_uy,a_uz,&
          a_vx,a_vy,a_vz,&
          a_wx,a_wy,a_wz,&
          mu,L,R,k) bind (C) 

       use iso_c_binding
       implicit none
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z
       real (c_double), value :: u_0
       real (c_double), value :: u_x
       real (c_double), value :: u_y
       real (c_double), value :: u_z
       real (c_double), value :: v_0
       real (c_double), value :: v_x
       real (c_double), value :: v_y
       real (c_double), value :: v_z
       real (c_double), value :: w_0
       real (c_double), value :: w_x
       real (c_double), value :: w_y
       real (c_double), value :: w_z
       real (c_double), value :: rho_0
       real (c_double), value :: rho_x
       real (c_double), value :: rho_y
       real (c_double), value :: rho_z
       real (c_double), value :: p_0
       real (c_double), value :: p_x
       real (c_double), value :: p_y
       real (c_double), value :: p_z
       real (c_double), value :: a_px
       real (c_double), value :: a_py
       real (c_double), value :: a_pz
       real (c_double), value :: a_rhox
       real (c_double), value :: a_rhoy
       real (c_double), value :: a_rhoz
       real (c_double), value :: a_ux
       real (c_double), value :: a_uy
       real (c_double), value :: a_uz
       real (c_double), value :: a_vx
       real (c_double), value :: a_vy
       real (c_double), value :: a_vz
       real (c_double), value :: a_wx
       real (c_double), value :: a_wy
       real (c_double), value :: a_wz
       real (c_double), value :: mu
       real (c_double), value :: L
       real (c_double), value :: R
       real (c_double), value :: k
       
     end function eval_3d_v_source

     real (c_double) function eval_3d_w_source(x,y,z,&
          u_0,u_x,u_y,u_z,&
          v_0,v_x,v_y,v_z,&
          w_0,w_x,w_y,w_z,&
          rho_0,rho_x,rho_y,rho_z,&
          p_0,p_x,p_y,p_z, &
          a_px,a_py,a_pz,&
          a_rhox,a_rhoy,a_rhoz,&
          a_ux,a_uy,a_uz,&
          a_vx,a_vy,a_vz,&
          a_wx,a_wy,a_wz,&
          mu,L,R,k) bind (C) 

       use iso_c_binding
       implicit none
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z
       real (c_double), value :: u_0
       real (c_double), value :: u_x
       real (c_double), value :: u_y
       real (c_double), value :: u_z
       real (c_double), value :: v_0
       real (c_double), value :: v_x
       real (c_double), value :: v_y
       real (c_double), value :: v_z
       real (c_double), value :: w_0
       real (c_double), value :: w_x
       real (c_double), value :: w_y
       real (c_double), value :: w_z
       real (c_double), value :: rho_0
       real (c_double), value :: rho_x
       real (c_double), value :: rho_y
       real (c_double), value :: rho_z
       real (c_double), value :: p_0
       real (c_double), value :: p_x
       real (c_double), value :: p_y
       real (c_double), value :: p_z
       real (c_double), value :: a_px
       real (c_double), value :: a_py
       real (c_double), value :: a_pz
       real (c_double), value :: a_rhox
       real (c_double), value :: a_rhoy
       real (c_double), value :: a_rhoz
       real (c_double), value :: a_ux
       real (c_double), value :: a_uy
       real (c_double), value :: a_uz
       real (c_double), value :: a_vx
       real (c_double), value :: a_vy
       real (c_double), value :: a_vz
       real (c_double), value :: a_wx
       real (c_double), value :: a_wy
       real (c_double), value :: a_wz
       real (c_double), value :: mu
       real (c_double), value :: L
       real (c_double), value :: R
       real (c_double), value :: k
       
     end function eval_3d_w_source

     real (c_double) function eval_3d_rho_source(x,y,z,&
          u_0,u_x,u_y,u_z,&
          v_0,v_x,v_y,v_z,&
          w_0,w_x,w_y,w_z,&
          rho_0,rho_x,rho_y,rho_z,&
          p_0,p_x,p_y,p_z, &
          a_px,a_py,a_pz,&
          a_rhox,a_rhoy,a_rhoz,&
          a_ux,a_uy,a_uz,&
          a_vx,a_vy,a_vz,&
          a_wx,a_wy,a_wz,&          
          mu,L,R,k) bind (C) 

       use iso_c_binding
       implicit none
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z
       real (c_double), value :: u_0
       real (c_double), value :: u_x
       real (c_double), value :: u_y
       real (c_double), value :: u_z
       real (c_double), value :: v_0
       real (c_double), value :: v_x
       real (c_double), value :: v_y
       real (c_double), value :: v_z
       real (c_double), value :: w_0
       real (c_double), value :: w_x
       real (c_double), value :: w_y
       real (c_double), value :: w_z
       real (c_double), value :: rho_0
       real (c_double), value :: rho_x
       real (c_double), value :: rho_y
       real (c_double), value :: rho_z
       real (c_double), value :: p_0
       real (c_double), value :: p_x
       real (c_double), value :: p_y
       real (c_double), value :: p_z
       real (c_double), value :: a_px
       real (c_double), value :: a_py
       real (c_double), value :: a_pz
       real (c_double), value :: a_rhox
       real (c_double), value :: a_rhoy
       real (c_double), value :: a_rhoz
       real (c_double), value :: a_ux
       real (c_double), value :: a_uy
       real (c_double), value :: a_uz
       real (c_double), value :: a_vx
       real (c_double), value :: a_vy
       real (c_double), value :: a_vz
       real (c_double), value :: a_wx
       real (c_double), value :: a_wy
       real (c_double), value :: a_wz
       real (c_double), value :: mu
       real (c_double), value :: L
       real (c_double), value :: R
       real (c_double), value :: k

     end function eval_3d_rho_source

     real (c_double) function eval_3d_e_source(&
          x,y,z,&
          u_0,u_x,u_y,u_z,&
          v_0,v_x,v_y,v_z,&
          w_0,w_x,w_y,w_z,&
          rho_0,rho_x,rho_y,rho_z,&
          p_0,p_x,p_y,p_z, &
          a_px,a_py,a_pz,&
          a_rhox,a_rhoy,a_rhoz,&
          a_ux,a_uy,a_uz,&
          a_vx,a_vy,a_vz,&
          a_wx,a_wy,a_wz,&          
          mu,Gamma,L,R,k) bind (C) 

       use iso_c_binding
       implicit none
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z
       real (c_double), value :: u_0
       real (c_double), value :: u_x
       real (c_double), value :: u_y
       real (c_double), value :: u_z
       real (c_double), value :: v_0
       real (c_double), value :: v_x
       real (c_double), value :: v_y
       real (c_double), value :: v_z
       real (c_double), value :: w_0
       real (c_double), value :: w_x
       real (c_double), value :: w_y
       real (c_double), value :: w_z
       real (c_double), value :: rho_0
       real (c_double), value :: rho_x
       real (c_double), value :: rho_y
       real (c_double), value :: rho_z
       real (c_double), value :: p_0
       real (c_double), value :: p_x
       real (c_double), value :: p_y
       real (c_double), value :: p_z
       real (c_double), value :: a_px
       real (c_double), value :: a_py
       real (c_double), value :: a_pz
       real (c_double), value :: a_rhox
       real (c_double), value :: a_rhoy
       real (c_double), value :: a_rhoz
       real (c_double), value :: a_ux
       real (c_double), value :: a_uy
       real (c_double), value :: a_uz
       real (c_double), value :: a_vx
       real (c_double), value :: a_vy
       real (c_double), value :: a_vz
       real (c_double), value :: a_wx
       real (c_double), value :: a_wy
       real (c_double), value :: a_wz
       real (c_double), value :: mu
       real (c_double), value :: Gamma
       real (c_double), value :: L
       real (c_double), value :: R
       real (c_double), value :: k

     end function eval_3d_e_source

     real (c_double) function eval_3d_p_an(x,y,z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,L) bind (C)
       use iso_c_binding
       implicit none

       real (c_double), value :: x,y,z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,L

     end function eval_3d_p_an
     
     real (c_double) function eval_3d_u_an(x,y,z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,L) bind (C)
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x,y,z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,L
       
     end function eval_3d_u_an

     real (c_double) function eval_3d_v_an(x,y,z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,L) bind (C)
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x,y,z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,L

     end function eval_3d_v_an

     real (c_double) function eval_3d_w_an(x,y,z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,L) bind (C)
       use iso_c_binding
       implicit none

       real (c_double), value :: x,y,z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,L
       
     end function eval_3d_w_an

     real (c_double) function eval_3d_rho_an(x,y,z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,L) bind (C)
       use iso_c_binding
       implicit none

       real (c_double), value :: x,y,z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,L
       
     end function eval_3d_rho_an

     !********************************************************
     !*
     !* 2D
     !*
     !********************************************************

     real (c_double) function eval_2d_u_source(x,y,&
          u_0,u_x,u_y,v_0,v_x,v_y,&
          rho_0,rho_x,rho_y,&
          p_0,p_x,p_y,a_px,a_py,&
          a_rhox,a_rhoy,&
          a_ux,a_uy,a_vx,a_vy,&
          mu,L,R,k) bind (C) 
       use iso_c_binding
       implicit none

       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: u_0
       real (c_double), value :: u_x
       real (c_double), value :: u_y
       real (c_double), value :: v_0
       real (c_double), value :: v_x
       real (c_double), value :: v_y
       real (c_double), value :: rho_0
       real (c_double), value :: rho_x
       real (c_double), value :: rho_y
       real (c_double), value :: p_0
       real (c_double), value :: p_x
       real (c_double), value :: p_y
       real (c_double), value :: a_px
       real (c_double), value :: a_py
       real (c_double), value :: a_rhox
       real (c_double), value :: a_rhoy
       real (c_double), value :: a_ux
       real (c_double), value :: a_uy
       real (c_double), value :: a_vx
       real (c_double), value :: a_vy
       real (c_double), value :: mu
       real (c_double), value :: L
       real (c_double), value :: R
       real (c_double), value :: k
              
     end function eval_2d_u_source

     real (c_double) function eval_2d_v_source(x,y,u_0,u_x,u_y,v_0,v_x,v_y, &
          rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py, &
          a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,mu,L,R,k) bind (C)
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: u_0
       real (c_double), value :: u_x
       real (c_double), value :: u_y
       real (c_double), value :: v_0
       real (c_double), value :: v_x
       real (c_double), value :: v_y
       real (c_double), value :: rho_0
       real (c_double), value :: rho_x
       real (c_double), value :: rho_y
       real (c_double), value :: p_0
       real (c_double), value :: p_x
       real (c_double), value :: p_y
       real (c_double), value :: a_px
       real (c_double), value :: a_py
       real (c_double), value :: a_rhox
       real (c_double), value :: a_rhoy
       real (c_double), value :: a_ux
       real (c_double), value :: a_uy
       real (c_double), value :: a_vx
       real (c_double), value :: a_vy
       real (c_double), value :: mu
       real (c_double), value :: L
       real (c_double), value :: R
       real (c_double), value :: k

     end function eval_2d_v_source

     real (c_double) function eval_2d_rho_source(x,y,u_0,u_x,u_y,v_0,v_x,v_y, &
          rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,mu,L,R,k) bind (C)
       use iso_c_binding
       implicit none

       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: u_0
       real (c_double), value :: u_x
       real (c_double), value :: u_y
       real (c_double), value :: v_0
       real (c_double), value :: v_x
       real (c_double), value :: v_y
       real (c_double), value :: rho_0
       real (c_double), value :: rho_x
       real (c_double), value :: rho_y
       real (c_double), value :: p_0
       real (c_double), value :: p_x
       real (c_double), value :: p_y
       real (c_double), value :: a_px
       real (c_double), value :: a_py
       real (c_double), value :: a_rhox
       real (c_double), value :: a_rhoy
       real (c_double), value :: a_ux
       real (c_double), value :: a_uy
       real (c_double), value :: a_vx
       real (c_double), value :: a_vy
       real (c_double), value :: mu
       real (c_double), value :: L
       real (c_double), value :: R
       real (c_double), value :: k

     end function eval_2d_rho_source

     real (c_double) function eval_2d_e_source(x,y,u_0,u_x,u_y,v_0,v_x,v_y, &
          rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy, &
          a_ux,a_uy,a_vx,a_vy,mu,L,R,k) bind (C)
       use iso_c_binding
       implicit none

       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: u_0
       real (c_double), value :: u_x
       real (c_double), value :: u_y
       real (c_double), value :: v_0
       real (c_double), value :: v_x
       real (c_double), value :: v_y
       real (c_double), value :: rho_0
       real (c_double), value :: rho_x
       real (c_double), value :: rho_y
       real (c_double), value :: p_0
       real (c_double), value :: p_x
       real (c_double), value :: p_y
       real (c_double), value :: a_px
       real (c_double), value :: a_py
       real (c_double), value :: a_rhox
       real (c_double), value :: a_rhoy
       real (c_double), value :: a_ux
       real (c_double), value :: a_uy
       real (c_double), value :: a_vx
       real (c_double), value :: a_vy
       real (c_double), value :: mu
       real (c_double), value :: L
       real (c_double), value :: R
       real (c_double), value :: k

     end function eval_2d_e_source

     real (c_double) function eval_2d_p_an(x,y,p_0,p_x,p_y,a_px,a_py,L) bind (C)
       use iso_c_binding
       implicit none

       real (c_double), value :: x,y,p_0,p_x,p_y,a_px,a_py,L

     end function eval_2d_p_an

     real (c_double) function eval_2d_u_an(x,y,u_0,u_x,u_y,a_ux,a_uy,L) bind (C)
       use iso_c_binding
       implicit none

       real (c_double), value :: x,y,u_0,u_x,u_y,a_ux,a_uy,L

     end function eval_2d_u_an

     real (c_double) function eval_2d_v_an(x,y,v_0,v_x,v_y,a_vx,a_vy,L) bind (C)
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x,y,v_0,v_x,v_y,a_vx,a_vy,L
     end function eval_2d_v_an

     real (c_double) function eval_2d_rho_an(x,y,rho_0,rho_x,rho_y,a_rhox,a_rhoy,L) bind (C)
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x,y,rho_0,rho_x,rho_y,a_rhox,a_rhoy,L
     end function eval_2d_rho_an

     !********************************************************
     !*
     !* 1D
     !*
     !********************************************************

     real (c_double) function eval_1d_u_source(x,u_0,u_x,rho_0,rho_x,p_0,p_x,a_px,a_rhox,a_ux,L) bind (C) 
       use iso_c_binding
       implicit none

       real (c_double), value :: x
       real (c_double), value :: u_0
       real (c_double), value :: u_x
       real (c_double), value :: rho_0
       real (c_double), value :: rho_x
       real (c_double), value :: p_0
       real (c_double), value :: p_x
       real (c_double), value :: a_px
       real (c_double), value :: a_rhox
       real (c_double), value :: a_ux
       real (c_double), value :: L
       
     end function eval_1d_u_source

     real (c_double) function eval_1d_rho_source(x,u_0,u_x,rho_0,rho_x,p_0,p_x,a_px,a_rhox,a_ux,L) bind (C) 
       use iso_c_binding
       implicit none

       real (c_double), value :: x
       real (c_double), value :: u_0
       real (c_double), value :: u_x
       real (c_double), value :: rho_0
       real (c_double), value :: rho_x
       real (c_double), value :: p_0
       real (c_double), value :: p_x
       real (c_double), value :: a_px
       real (c_double), value :: a_rhox
       real (c_double), value :: a_ux
       real (c_double), value :: L

     end function eval_1d_rho_source

     real (c_double) function eval_1d_e_source(x,u_0,u_x,rho_0,rho_x,p_0,p_x,a_px,a_rhox,a_ux,Gamma,mu,L) bind (C) 
       use iso_c_binding
       implicit none

       real (c_double), value :: x
       real (c_double), value :: u_0
       real (c_double), value :: u_x
       real (c_double), value :: rho_0
       real (c_double), value :: rho_x
       real (c_double), value :: p_0
       real (c_double), value :: p_x
       real (c_double), value :: a_px
       real (c_double), value :: a_rhox
       real (c_double), value :: a_ux
       real (c_double), value :: mu
       real (c_double), value :: Gamma
       real (c_double), value :: L

     end function eval_1d_e_source

     real (c_double) function eval_1d_p_an(x,p_0,p_x,a_px,L) bind (C)
       use iso_c_binding
       implicit none

       real (c_double), value :: x,p_0,p_x,a_px,L

     end function eval_1d_p_an

     real (c_double) function eval_1d_u_an(x,p_0,p_x,a_px,L) bind (C)
       use iso_c_binding
       implicit none

       real (c_double), value :: x,p_0,p_x,a_px,L

     end function eval_1d_u_an
     
     real (c_double) function eval_1d_v_an(x,p_0,p_x,a_px,L) bind (C)
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x,p_0,p_x,a_px,L
     end function eval_1d_v_an

     real (c_double) function eval_1d_rho_an(x,p_0,p_x,a_px,L) bind (C)
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x,p_0,p_x,a_px,L
     end function eval_1d_rho_an

  end interface
       
end module cns_source_interface

