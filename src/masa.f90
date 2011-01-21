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
! masa.f90: Fortran module interface definition for MASA routines
!
! $Id$
! --------------------------------------------------------------------------
! --------------------------------------------------------------------------

module masa

  ! ---------------------------------
  ! MMS Init/Selection Routines
  ! ---------------------------------

  interface
     subroutine masa_init_passthrough(user_tag,desired_mms_function) bind (C,name='cmasa_init')
       use iso_c_binding
       implicit none

       character(c_char), intent(in) :: user_tag(*)
       character(c_char), intent(in) :: desired_mms_function(*)

     end subroutine masa_init_passthrough
  end interface

  interface
     subroutine masa_list_mms() bind (C,name='cmasa_list_mms')
       use iso_c_binding
       implicit none

     end subroutine masa_list_mms
  end interface

  interface
     subroutine masa_select_mms_passthrough(desired_mms_function) bind (C,name='cmasa_select_mms')
       use iso_c_binding
       implicit none

       character(c_char), intent(in) :: desired_mms_function(*)

     end subroutine masa_select_mms_passthrough
  end interface

  interface
     subroutine masa_sanity_check() bind (C,name='cmasa_sanity_check')
       use iso_c_binding
       implicit none

     end subroutine masa_sanity_check
  end interface

  ! ---------------------------------
  ! MMS Parameter Routines
  ! ---------------------------------

  interface
     subroutine masa_init_param() bind (C,name='cmasa_init_param')
       use iso_c_binding
       implicit none

     end subroutine masa_init_param
  end interface

  interface
     subroutine masa_display_param() bind (C,name='cmasa_display_param')
       use iso_c_binding
       implicit none

     end subroutine masa_display_param
  end interface

  interface
     real (c_double) function masa_get_param_passthrough(param_name) bind (C,name='cmasa_get_param')
       use iso_c_binding
       implicit none

       character(c_char), intent(in) :: param_name(*)

     end function masa_get_param_passthrough
  end interface

  interface
     subroutine masa_set_param_passthrough(param_name,value) bind (C,name='cmasa_set_param')
       use iso_c_binding
       implicit none

       character(c_char), intent(in) :: param_name(*)
       real (c_double), value        :: value

     end subroutine masa_set_param_passthrough
  end interface

  ! ---------------------------------
  ! MMS source term interfaces -- 1d
  ! ---------------------------------

  interface 
     real (c_double) function masa_eval_1d_t_source(value) bind (C,name='cmasa_eval_1d_t_source')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_t_source
  end interface

  interface 
     real (c_double) function masa_eval_1d_u_source(value) bind (C,name='cmasa_eval_1d_u_source')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_u_source
  end interface

  interface 
     real (c_double) function masa_eval_1d_e_source(value) bind (C,name='cmasa_eval_1d_e_source')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_e_source
  end interface

  interface 
     real (c_double) function masa_eval_1d_rho_source(value) bind (C,name='cmasa_eval_1d_rho_source')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_rho_source
  end interface  

  ! ---------------------------------
  ! MMS source term interfaces -- 2d
  ! ---------------------------------

  interface 
     real (c_double) function masa_eval_2d_t_source(value,value2) bind (C,name='cmasa_eval_2d_t_source')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_t_source
  end interface

  interface 
     real (c_double) function masa_eval_2d_u_source(value,value2) bind (C,name='cmasa_eval_2d_u_source')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2

     end function masa_eval_2d_u_source
  end interface

  interface 
     real (c_double) function masa_eval_2d_v_source(value,value2) bind (C,name='cmasa_eval_2d_v_source')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_v_source
  end interface

  interface 
     real (c_double) function masa_eval_2d_e_source(value,value2) bind (C,name='cmasa_eval_2d_e_source')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_e_source
  end interface

  interface 
     real (c_double) function masa_eval_2d_rho_source(value,value2) bind (C,name='cmasa_eval_2d_rho_source')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_rho_source
  end interface  

  ! ---------------------------------
  ! MMS source term interfaces -- 3d
  ! ---------------------------------

  interface 
     real (c_double) function masa_eval_3d_t_source(value,value2,value3) bind (C,name='cmasa_eval_3d_t_source')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3
       
     end function masa_eval_3d_t_source
  end interface

  interface 
     real (c_double) function masa_eval_3d_u_source(value,value2,value3) bind (C,name='cmasa_eval_3d_u_source')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3
       
     end function masa_eval_3d_u_source
  end interface

  interface 
     real (c_double) function masa_eval_3d_v_source(value,value2,value3) bind (C,name='cmasa_eval_3d_v_source')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       

     end function masa_eval_3d_v_source
  end interface

  interface 
     real (c_double) function masa_eval_3d_w_source(value,value2,value3) bind (C,name='cmasa_eval_3d_w_source')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       

     end function masa_eval_3d_w_source
  end interface

  interface 
     real (c_double) function masa_eval_3d_e_source(value,value2,value3) bind (C,name='cmasa_eval_3d_e_source')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3
       
     end function masa_eval_3d_e_source
  end interface

  interface 
     real (c_double) function masa_eval_3d_rho_source(value,value2,value3) bind (C,name='cmasa_eval_3d_rho_source')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       

     end function masa_eval_3d_rho_source
  end interface  

  ! ---------------------------------
  ! MMS analytical term interfaces -- 1d
  ! ---------------------------------

  interface 
     real (c_double) function masa_eval_1d_t_an(value) bind (C,name='cmasa_eval_1d_t_an')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_t_an
  end interface

  interface 
     real (c_double) function masa_eval_1d_u_an(value) bind (C,name='cmasa_eval_1d_u_an')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_u_an
  end interface

  interface 
     real (c_double) function masa_eval_1d_p_an(value) bind (C,name='cmasa_eval_1d_p_an')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_p_an
  end interface

  interface 
     real (c_double) function masa_eval_1d_rho_an(value) bind (C,name='cmasa_eval_1d_rho_an')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_rho_an
  end interface

  ! ---------------------------------
  ! MMS analytical term interfaces -- 2d
  ! ---------------------------------

  interface 
     real (c_double) function masa_eval_2d_t_an(value,value2) bind (C,name='cmasa_eval_2d_t_an')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_t_an
  end interface

  interface 
     real (c_double) function masa_eval_2d_u_an(value,value2) bind (C,name='cmasa_eval_2d_u_an')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_u_an
  end interface

  interface 
     real (c_double) function masa_eval_2d_v_an(value,value2) bind (C,name='cmasa_eval_2d_v_an')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_v_an
  end interface

  interface 
     real (c_double) function masa_eval_2d_p_an(value,value2) bind (C,name='cmasa_eval_2d_p_an')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_p_an
  end interface

  interface 
     real (c_double) function masa_eval_2d_rho_an(value,value2) bind (C,name='cmasa_eval_2d_rho_an')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_rho_an
  end interface

  ! ---------------------------------
  ! MMS analytical term interfaces -- 3d
  ! ---------------------------------

  interface 
     real (c_double) function masa_eval_3d_t_an(value,value2,value3) bind (C,name='cmasa_eval_3d_t_an')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       

     end function masa_eval_3d_t_an
  end interface

  interface 
     real (c_double) function masa_eval_3d_u_an(value,value2,value3) bind (C,name='cmasa_eval_3d_u_an')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       
       
     end function masa_eval_3d_u_an
  end interface

  interface 
     real (c_double) function masa_eval_3d_v_an(value,value2,value3) bind (C,name='cmasa_eval_3d_v_an')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       
       
     end function masa_eval_3d_v_an
  end interface

  interface 
     real (c_double) function masa_eval_3d_w_an(value,value2,value3) bind (C,name='cmasa_eval_3d_w_an')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       
       
     end function masa_eval_3d_w_an
  end interface

  interface 
     real (c_double) function masa_eval_3d_p_an(value,value2,value3) bind (C,name='cmasa_eval_3d_p_an')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       
       
     end function masa_eval_3d_p_an
  end interface

  interface 
     real (c_double) function masa_eval_3d_rho_an(value,value2,value3) bind (C,name='cmasa_eval_3d_rho_an')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       
       
     end function masa_eval_3d_rho_an
  end interface

  ! ---------------------------------
  ! MMS gradient term interfaces -- 1d
  ! ---------------------------------

  interface 
     real (c_double) function masa_eval_1d_grad_u(value) bind (C,name='cmasa_eval_1d_grad_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_grad_u
  end interface

  interface 
     real (c_double) function masa_eval_1d_grad_p(value) bind (C,name='cmasa_eval_1d_grad_p')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_grad_p
  end interface

  interface 
     real (c_double) function masa_eval_1d_grad_rho(value) bind (C,name='cmasa_eval_1d_grad_rho')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_grad_rho
  end interface

  ! ---------------------------------
  ! MMS gradient term interfaces -- 2d
  ! ---------------------------------

  interface 
     real (c_double) function masa_eval_2d_grad_u(value,value2,it) bind (C,name='cmasa_eval_2d_grad_u')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: value
       real    (c_double), value :: value2
       integer (c_int),    value :: it
       
     end function masa_eval_2d_grad_u
  end interface

  interface 
     real (c_double) function masa_eval_2d_grad_v(value,value2,it) bind (C,name='cmasa_eval_2d_grad_v')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: value
       real    (c_double), value :: value2
       integer (c_int),    value :: it
       
     end function masa_eval_2d_grad_v
  end interface

  interface 
     real (c_double) function masa_eval_2d_grad_w(value,value2,it) bind (C,name='cmasa_eval_2d_grad_w')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: value
       real    (c_double), value :: value2
       integer (c_int),    value :: it
       
     end function masa_eval_2d_grad_w
  end interface

  interface 
     real (c_double) function masa_eval_2d_grad_p(value,value2,it) bind (C,name='cmasa_eval_2d_grad_p')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: value
       real    (c_double), value :: value2
       integer (c_int),    value :: it
       
     end function masa_eval_2d_grad_p
  end interface

  interface 
     real (c_double) function masa_eval_2d_grad_rho(value,value2,it) bind (C,name='cmasa_eval_2d_grad_rho')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: value
       real    (c_double), value :: value2
       integer (c_int),    value :: it
       
     end function masa_eval_2d_grad_rho
  end interface

  ! ---------------------------------
  ! MMS gradient term interfaces -- 3d
  ! ---------------------------------

  interface 
     real (c_double) function masa_eval_3d_grad_u(value,value2,value3,it) bind (C,name='cmasa_eval_3d_grad_u')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: value
       real    (c_double), value :: value2
       real    (c_double), value :: value3
       integer (c_int),    value :: it
       
     end function masa_eval_3d_grad_u
  end interface

  interface 
     real (c_double) function masa_eval_3d_grad_v(value,value2,value3,it) bind (C,name='cmasa_eval_3d_grad_v')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: value
       real    (c_double), value :: value2
       real    (c_double), value :: value3
       integer (c_int),    value :: it
       
     end function masa_eval_3d_grad_v
  end interface

  interface 
     real (c_double) function masa_eval_3d_grad_w(value,value2,value3,it) bind (C,name='cmasa_eval_3d_grad_w')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: value
       real    (c_double), value :: value2
       real    (c_double), value :: value3
       integer (c_int),    value :: it
       
     end function masa_eval_3d_grad_w
  end interface

  interface 
     real (c_double) function masa_eval_3d_grad_p(value,value2,value3,it) bind (C,name='cmasa_eval_3d_grad_p')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: value
       real    (c_double), value :: value2
       real    (c_double), value :: value3
       integer (c_int),    value :: it
       
     end function masa_eval_3d_grad_p
  end interface

  interface 
     real (c_double) function masa_eval_3d_grad_rho(value,value2,value3,it) bind (C,name='cmasa_eval_3d_grad_rho')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: value
       real    (c_double), value :: value2
       real    (c_double), value :: value3
       integer (c_int),    value :: it
       
     end function masa_eval_3d_grad_rho
  end interface
  
contains 
  
  ! ----------------------------------------------------------------
  ! Wrapper routines for functions which include character
  ! strings; the wrapers insert necessary null terminators for 
  ! subsequent use with C/C++
  ! ----------------------------------------------------------------
  
  subroutine masa_init(user_tag,desired_mms_function)
    use iso_c_binding
    implicit none

    character(len=*) :: user_tag
    character(len=*) :: desired_mms_function

    call masa_init_passthrough(user_tag//C_NULL_CHAR,desired_mms_function//C_NULL_CHAR)
    return
  end subroutine masa_init

  subroutine masa_select_mms(desired_mms_function)
    use iso_c_binding
    implicit none

    character(len=*) :: desired_mms_function

    call masa_select_mms_passthrough(desired_mms_function//C_NULL_CHAR)
    return
  end subroutine masa_select_mms

  real (c_double) function masa_get_param(param_name)
    use iso_c_binding
    implicit none

    character(len=*) :: param_name

    masa_get_param =  masa_get_param_passthrough(param_name//C_NULL_CHAR)

  end function masa_get_param

  subroutine masa_set_param(param_name,value)
    use iso_c_binding
    implicit none

    character(len=*), intent(in)        :: param_name
    real (c_double), intent(in), value  :: value

    call masa_set_param_passthrough(param_name//C_NULL_CHAR,value)

  end subroutine masa_set_param

end module masa
