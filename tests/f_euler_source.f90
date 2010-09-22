!! -*-f90-*-
!!-----------------------------------------------------------------------bl-
!!-----------------------------------------------------------------------el-
!! $Author: nick $
!! $Id$
!!--------------------------------------------------------------------------
!!--------------------------------------------------------------------------

module euler_source_interface

  interface
     real (c_double) function eval_2d_u_source(x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,&
          p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,L) bind (C) 
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
       real (c_double), value :: L
       
     end function eval_2d_u_source

     real (c_double) function eval_2d_v_source(x,y,u_0,u_x,u_y,v_0,v_x,v_y, &
          rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py, &
          a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,L) bind (C)
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
       real (c_double), value :: L
     end function eval_2d_v_source

     real (c_double) function eval_2d_rho_source(x,y,u_0,u_x,u_y,v_0,v_x,v_y, &
          rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py, &
          a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,L) bind (C)
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
       real (c_double), value :: L
     end function eval_2d_rho_source

     real (c_double) function eval_2d_e_source(x,y,u_0,u_x,u_y,v_0,v_x,v_y, &
          rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py, &
          a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,L) bind (C)
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
       real (c_double), value :: L
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

  end interface
       

end module euler_source_interface

