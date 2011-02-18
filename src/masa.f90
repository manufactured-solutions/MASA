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
  use iso_c_binding
  implicit none
  type (c_funptr),bind(c)          :: funp
    
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
     subroutine masa_purge_default_param() bind (C,name='cmasa_purge_default_param')
       use iso_c_binding
       implicit none
       
     end subroutine masa_purge_default_param
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
     real (c_double) function masa_eval_1d_source_t(value) bind (C,name='cmasa_eval_1d_source_t')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_source_t
  end interface

  interface 
     real (c_double) function masa_eval_1d_source_u(value) bind (C,name='cmasa_eval_1d_source_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_source_u
  end interface

  interface 
     real (c_double) function masa_eval_1d_source_e(value) bind (C,name='cmasa_eval_1d_source_e')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_source_e
  end interface

  interface 
     real (c_double) function masa_eval_1d_source_rho(value) bind (C,name='cmasa_eval_1d_source_rho')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_source_rho
  end interface  

  interface 
     real (c_double) function masa_eval_1d_source_rho_u(value) bind (C,name='cmasa_eval_1d_source_rho_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_source_rho_u
  end interface  

  interface 
     real (c_double) function masa_eval_1d_source_rho_e(value) bind (C,name='cmasa_eval_1d_source_rho_e')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_source_rho_e
  end interface

  interface 
     real (c_double) function masa_eval_1d_source_rho_N_passthrough(value,fun) bind (C,name='cmasa_eval_1d_source_rho_N')
       use iso_c_binding
       implicit none
       
       real (c_double), value      :: value
       type (c_funptr), intent(in) :: fun
       
     end function masa_eval_1d_source_rho_N_passthrough
  end interface  


  interface 
     real (c_double) function masa_eval_1d_source_rho_N2_passthrough(value,fun) bind (C,name='cmasa_eval_1d_source_rho_N2')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       type (c_funptr), intent(in) :: fun

     end function masa_eval_1d_source_rho_N2_passthrough
  end interface  

  ! ---------------------------------
  ! MMS source term interfaces -- 2d
  ! ---------------------------------

  interface 
     real (c_double) function masa_eval_2d_source_t(value,value2) bind (C,name='cmasa_eval_2d_source_t')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_source_t
  end interface

  interface 
     real (c_double) function masa_eval_2d_source_u(value,value2) bind (C,name='cmasa_eval_2d_source_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2

     end function masa_eval_2d_source_u
  end interface

  interface 
     real (c_double) function masa_eval_2d_source_v(value,value2) bind (C,name='cmasa_eval_2d_source_v')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_source_v
  end interface

  interface 
     real (c_double) function masa_eval_2d_source_e(value,value2) bind (C,name='cmasa_eval_2d_source_e')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_source_e
  end interface

  interface 
     real (c_double) function masa_eval_2d_source_rho(value,value2) bind (C,name='cmasa_eval_2d_source_rho')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_source_rho
  end interface  

  interface 
     real (c_double) function masa_eval_2d_source_rho_u(value,value2) bind (C,name='cmasa_eval_2d_source_rho_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_source_rho_u
  end interface  

  interface 
     real (c_double) function masa_eval_2d_source_rho_v(value,value2) bind (C,name='cmasa_eval_2d_source_rho_v')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_source_rho_v
  end interface
  
  interface 
     real (c_double) function masa_eval_2d_source_rho_w(value,value2) bind (C,name='cmasa_eval_2d_source_rho_w')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_source_rho_w
  end interface  

  interface 
     real (c_double) function masa_eval_2d_source_rho_e(value,value2) bind (C,name='cmasa_eval_2d_source_rho_e')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_source_rho_e
  end interface  

  ! ---------------------------------
  ! MMS source term interfaces -- 3d
  ! ---------------------------------

  interface 
     real (c_double) function masa_eval_3d_source_t(value,value2,value3) bind (C,name='cmasa_eval_3d_source_t')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3
       
     end function masa_eval_3d_source_t
  end interface

  interface 
     real (c_double) function masa_eval_3d_source_u(value,value2,value3) bind (C,name='cmasa_eval_3d_source_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3
       
     end function masa_eval_3d_source_u
  end interface

  interface 
     real (c_double) function masa_eval_3d_source_v(value,value2,value3) bind (C,name='cmasa_eval_3d_source_v')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       

     end function masa_eval_3d_source_v
  end interface

  interface 
     real (c_double) function masa_eval_3d_source_w(value,value2,value3) bind (C,name='cmasa_eval_3d_source_w')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       

     end function masa_eval_3d_source_w
  end interface

  interface 
     real (c_double) function masa_eval_3d_source_e(value,value2,value3) bind (C,name='cmasa_eval_3d_source_e')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3
       
     end function masa_eval_3d_source_e
  end interface

  interface 
     real (c_double) function masa_eval_3d_source_rho(value,value2,value3) bind (C,name='cmasa_eval_3d_source_rho')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       

     end function masa_eval_3d_source_rho
  end interface  

  interface 
     real (c_double) function masa_eval_3d_source_rho_u(value,value2,value3) bind (C,name='cmasa_eval_3d_source_rho_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       

     end function masa_eval_3d_source_rho_u
  end interface  

  interface 
     real (c_double) function masa_eval_3d_source_rho_v(value,value2,value3) bind (C,name='cmasa_eval_3d_source_rho_v')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       

     end function masa_eval_3d_source_rho_v
  end interface  

  interface 
     real (c_double) function masa_eval_3d_source_rho_w(value,value2,value3) bind (C,name='cmasa_eval_3d_source_rho_w')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       

     end function masa_eval_3d_source_rho_w
  end interface  

  interface 
     real (c_double) function masa_eval_3d_source_rho_e(value,value2,value3) bind (C,name='cmasa_eval_3d_source_rho_e')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       
       
     end function masa_eval_3d_source_rho_e
  end interface  

  ! ---------------------------------
  ! MMS analytical term interfaces -- 1d
  ! ---------------------------------

  interface 
     real (c_double) function masa_eval_1d_exact_t(value) bind (C,name='cmasa_eval_1d_exact_t')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_exact_t
  end interface

  interface 
     real (c_double) function masa_eval_1d_exact_u(value) bind (C,name='cmasa_eval_1d_exact_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_exact_u
  end interface

  interface 
     real (c_double) function masa_eval_1d_exact_p(value) bind (C,name='cmasa_eval_1d_exact_p')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_exact_p
  end interface

  interface 
     real (c_double) function masa_eval_1d_exact_rho(value) bind (C,name='cmasa_eval_1d_exact_rho')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_exact_rho
  end interface

  interface 
     real (c_double) function masa_eval_1d_exact_rho_N(value) bind (C,name='cmasa_eval_1d_exact_rho_N')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_exact_rho_N
  end interface

  interface 
     real (c_double) function masa_eval_1d_exact_rho_N2(value) bind (C,name='cmasa_eval_1d_exact_rho_N2')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_1d_exact_rho_N2
  end interface

  ! ---------------------------------
  ! MMS analytical term interfaces -- 2d
  ! ---------------------------------

  interface 
     real (c_double) function masa_eval_2d_exact_t(value,value2) bind (C,name='cmasa_eval_2d_exact_t')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_exact_t
  end interface

  interface 
     real (c_double) function masa_eval_2d_exact_u(value,value2) bind (C,name='cmasa_eval_2d_exact_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_exact_u
  end interface

  interface 
     real (c_double) function masa_eval_2d_exact_v(value,value2) bind (C,name='cmasa_eval_2d_exact_v')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_exact_v
  end interface

  interface 
     real (c_double) function masa_eval_2d_exact_p(value,value2) bind (C,name='cmasa_eval_2d_exact_p')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_exact_p
  end interface

  interface 
     real (c_double) function masa_eval_2d_exact_rho(value,value2) bind (C,name='cmasa_eval_2d_exact_rho')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       
     end function masa_eval_2d_exact_rho
  end interface

  ! ---------------------------------
  ! MMS analytical term interfaces -- 3d
  ! ---------------------------------

  interface 
     real (c_double) function masa_eval_3d_exact_t(value,value2,value3) bind (C,name='cmasa_eval_3d_exact_t')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       

     end function masa_eval_3d_exact_t
  end interface

  interface 
     real (c_double) function masa_eval_3d_exact_u(value,value2,value3) bind (C,name='cmasa_eval_3d_exact_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       
       
     end function masa_eval_3d_exact_u
  end interface

  interface 
     real (c_double) function masa_eval_3d_exact_v(value,value2,value3) bind (C,name='cmasa_eval_3d_exact_v')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       
       
     end function masa_eval_3d_exact_v
  end interface

  interface 
     real (c_double) function masa_eval_3d_exact_w(value,value2,value3) bind (C,name='cmasa_eval_3d_exact_w')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       
       
     end function masa_eval_3d_exact_w
  end interface

  interface 
     real (c_double) function masa_eval_3d_exact_p(value,value2,value3) bind (C,name='cmasa_eval_3d_exact_p')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       
       
     end function masa_eval_3d_exact_p
  end interface

  interface 
     real (c_double) function masa_eval_3d_exact_rho(value,value2,value3) bind (C,name='cmasa_eval_3d_exact_rho')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       real (c_double), value :: value2
       real (c_double), value :: value3       
       
     end function masa_eval_3d_exact_rho
  end interface

  ! ---------------------------------
  ! MMS gradient term interfaces -- 1d
  ! ---------------------------------

  interface 
     real (c_double) function masa_eval_grad_1d_u(value) bind (C,name='cmasa_eval_grad_1d_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_grad_1d_u
  end interface

  interface 
     real (c_double) function masa_eval_grad_1d_p(value) bind (C,name='cmasa_eval_grad_1d_p')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_grad_1d_p
  end interface

  interface 
     real (c_double) function masa_eval_grad_1d_rho(value) bind (C,name='cmasa_eval_grad_1d_rho')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: value
       
     end function masa_eval_grad_1d_rho
  end interface

  ! ---------------------------------
  ! MMS gradient term interfaces -- 2d
  ! ---------------------------------

  interface 
     real (c_double) function masa_eval_grad_2d_u(value,value2,it) bind (C,name='cmasa_eval_grad_2d_u')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: value
       real    (c_double), value :: value2
       integer (c_int),    value :: it
       
     end function masa_eval_grad_2d_u
  end interface

  interface 
     real (c_double) function masa_eval_grad_2d_v(value,value2,it) bind (C,name='cmasa_eval_grad_2d_v')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: value
       real    (c_double), value :: value2
       integer (c_int),    value :: it
       
     end function masa_eval_grad_2d_v
  end interface

  interface 
     real (c_double) function masa_eval_grad_2d_w(value,value2,it) bind (C,name='cmasa_eval_grad_2d_w')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: value
       real    (c_double), value :: value2
       integer (c_int),    value :: it
       
     end function masa_eval_grad_2d_w
  end interface

  interface 
     real (c_double) function masa_eval_grad_2d_p(value,value2,it) bind (C,name='cmasa_eval_grad_2d_p')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: value
       real    (c_double), value :: value2
       integer (c_int),    value :: it
       
     end function masa_eval_grad_2d_p
  end interface

  interface 
     real (c_double) function masa_eval_grad_2d_rho(value,value2,it) bind (C,name='cmasa_eval_grad_2d_rho')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: value
       real    (c_double), value :: value2
       integer (c_int),    value :: it
       
     end function masa_eval_grad_2d_rho
  end interface

  ! ---------------------------------
  ! MMS gradient term interfaces -- 3d
  ! ---------------------------------

  interface 
     real (c_double) function masa_eval_grad_3d_u(value,value2,value3,it) bind (C,name='cmasa_eval_grad_3d_u')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: value
       real    (c_double), value :: value2
       real    (c_double), value :: value3
       integer (c_int),    value :: it
       
     end function masa_eval_grad_3d_u
  end interface

  interface 
     real (c_double) function masa_eval_grad_3d_v(value,value2,value3,it) bind (C,name='cmasa_eval_grad_3d_v')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: value
       real    (c_double), value :: value2
       real    (c_double), value :: value3
       integer (c_int),    value :: it
       
     end function masa_eval_grad_3d_v
  end interface

  interface 
     real (c_double) function masa_eval_grad_3d_w(value,value2,value3,it) bind (C,name='cmasa_eval_grad_3d_w')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: value
       real    (c_double), value :: value2
       real    (c_double), value :: value3
       integer (c_int),    value :: it
       
     end function masa_eval_grad_3d_w
  end interface

  interface 
     real (c_double) function masa_eval_grad_3d_p(value,value2,value3,it) bind (C,name='cmasa_eval_grad_3d_p')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: value
       real    (c_double), value :: value2
       real    (c_double), value :: value3
       integer (c_int),    value :: it
       
     end function masa_eval_grad_3d_p
  end interface

  interface 
     real (c_double) function masa_eval_grad_3d_rho(value,value2,value3,it) bind (C,name='cmasa_eval_grad_3d_rho')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: value
       real    (c_double), value :: value2
       real    (c_double), value :: value3
       integer (c_int),    value :: it
       
     end function masa_eval_grad_3d_rho
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

  real (c_double) function masa_eval_1d_source_rho_N(value,funct)
    use iso_c_binding
    implicit none
    
    real(c_double), intent(in)       :: value
    real(c_double), external,bind(c) :: funct
    print *, "here !!"
    
    funp = c_funloc(funct)
    masa_eval_1d_source_rho_N = masa_eval_1d_source_rho_N_passthrough(value,funp)

  end function masa_eval_1d_source_rho_N
  
  real (c_double) function masa_eval_1d_source_rho_N2(value,functionguy)
    use iso_c_binding
    implicit none
    
    real(c_double), intent(in)       :: value
    real(c_double), external,bind(c) :: functionguy
    funp = c_funloc(functionguy)
    masa_eval_1d_source_rho_N2 = masa_eval_1d_source_rho_N2_passthrough(value,funp)
  
  end function masa_eval_1d_source_rho_N2

end module masa
