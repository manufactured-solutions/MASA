!! -*-f90-*-
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

  !> \file
  !! MASA Fortran Interface

  ! -------------------------------------
  !! \name MMS Init/Selection Routines
  ! -------------------------------------

  interface
     !> Initalizes a masa manufactured solution class.
     !! @param user_tag The first character string is a handle 
     !! for the newly initalized class. (e.g. "nick's mms")
     !!
     !! @param mms_id The second character string is the unique masa 
     !! identifier string for a particular masa class. (e.g. "euler_1d"). 
     !! See \subpage mms_avail "Available Manufactured Solutions"
     !! for information on all available solution strings.
     !!
     subroutine masa_init_passthrough(user_tag,mms_id) bind (C,name='masa_init')
       use iso_c_binding
       implicit none

       character(c_char), intent(in) :: user_tag(*)
       character(c_char), intent(in) :: mms_id(*)

     end subroutine masa_init_passthrough
  end interface

  interface
     !> Display (to stdout) the number of user initalized solutions. 
     !!
     !! This will then display each solutions unique handle and 
     !! full manufactured class name.
     !!
     subroutine masa_list_mms() bind (C,name='masa_list_mms')
       use iso_c_binding
       implicit none

     end subroutine masa_list_mms
  end interface

  interface     
     !> This function sets all masa parameters of the currently 
     !! initialized solution to uninitalized.
     !!
     subroutine masa_purge_default_param() bind (C,name='masa_purge_default_param')
       use iso_c_binding
       implicit none
       
     end subroutine masa_purge_default_param
  end interface

  interface
     !> Selects an already initalized manufactured solution class.
     !!
     !! Thus, if the user had created two manufactured 
     !! classes (say, nick and bob) using masa_init, 
     !! he could switch between them by passing the 
     !! handle to this routine. ex. masa_select_mms("nick")
     !!
     !! @param desired_mms_function Unique manufactured class handle string.
     !! 
     subroutine masa_select_mms_passthrough(desired_mms_function) bind (C,name='masa_select_mms')
       use iso_c_binding
       implicit none

       character(c_char), intent(in) :: desired_mms_function(*)

     end subroutine masa_select_mms_passthrough
  end interface

  interface
     !> Checks that all parameters for the currently selected 
     !! manufactured class have been initalized to some value.
     !!
     subroutine masa_sanity_check() bind (C,name='masa_sanity_check')
       use iso_c_binding
       implicit none

     end subroutine masa_sanity_check
  end interface

  ! ---------------------------------
  !! /name MMS Parameter Routines
  ! ---------------------------------

  interface
     !> Subroutine that will initalize
     !! all the registered variables to selected defaults
     !! for the currently selected manufactured solution class. 
     !!
     subroutine masa_init_param() bind (C,name='masa_init_param')
       use iso_c_binding
       implicit none

     end subroutine masa_init_param
  end interface

  interface
     !> Output the currently selected manufactured
     !! solution class' parameter names and values to standard output
     !!
     subroutine masa_display_param() bind (C,name='masa_display_param')
       use iso_c_binding
       implicit none
       
     end subroutine masa_display_param
  end interface

  interface
     
     !> Will return a particular registered variables inside 
     !! the currently selected manufactured solution class. 
     !!
     !! @param param_name a character string for
     !! the particular variable value to be returned.
     !!
     !! @return a double of the currently held value 
     !! of the specified variable. 
     !!
     real (c_double) function masa_get_param_passthrough(param_name) bind (C,name='masa_get_param')
       use iso_c_binding
       implicit none

       character(c_char), intent(in) :: param_name(*)

     end function masa_get_param_passthrough
  end interface

  interface
     !> Set a particular registered variables inside 
     !! the currently selected manufactured solution class. 
     !!
     !! @param param_name A character string for the 
     !! particular variable to be set. 
     !!
     !! @param Real(8) number to use as the new 
     !! value of the variable.
     !!
     subroutine masa_set_param_passthrough(param_name,value) bind (C,name='masa_set_param')
       use iso_c_binding
       implicit none

       character(c_char), intent(in) :: param_name(*)
       real (c_double), value        :: value

     end subroutine masa_set_param_passthrough
  end interface

  ! ---------------------------------
  ! MMS Vector/Array Routines
  ! ---------------------------------

  interface
     !> Output the currently selected manufactured
     !! solution class' vector names and lengths to standard output.
     !!
     subroutine masa_display_array() bind (C,name='masa_display_array')
       use iso_c_binding
       implicit none

     end subroutine masa_display_array
  end interface

  interface
     !> Return an array from the currently selected manufactured 
     !! solution class. 
     !!
     !! @param array_name Name of the desired array.
     !!
     !! @param it Integer specifying the length of the array.
     !!
     !! @param array The array.
     !!
     subroutine masa_get_array_passthrough(array_name,it,array) bind (C,name='masa_get_array')
       use iso_c_binding
       implicit none

       character(c_char), intent(in)                :: array_name(*)
       integer  (c_int)                             :: it      ! pass-by-ref is intentional
       real     (c_double),dimension(*),intent(out) :: array

     end subroutine masa_get_array_passthrough
  end interface

!  interface
!     subroutine masa_set_array_passthrough(param_name,value) bind (C,name='masa_set_array')
!       use iso_c_binding
!       implicit none
!       
!       character(c_char), intent(in) :: param_name(*)
!       real (c_double), value        :: value
!       
!     end subroutine masa_set_array_passthrough
!  end interface

  ! ---------------------------------
  ! MMS source term interfaces -- 1d
  ! ---------------------------------

  interface 
     !> Evaluates the one dimensional source term of the temperature.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_1d_source_t(x) bind (C,name='masa_eval_1d_source_t')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       
     end function masa_eval_1d_source_t
  end interface

  interface 
     !> Evaluates the one dimensional source term of the velocity.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_1d_source_u(x) bind (C,name='masa_eval_1d_source_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       
     end function masa_eval_1d_source_u
  end interface

  interface 
     !> Evaluates the one dimensional source term of the energy.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_1d_source_e(x) bind (C,name='masa_eval_1d_source_e')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       
     end function masa_eval_1d_source_e
  end interface

  interface 
     !> Evaluates the one dimensional source term of the density.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_1d_source_rho(x) bind (C,name='masa_eval_1d_source_rho')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       
     end function masa_eval_1d_source_rho
  end interface  

  interface 
     !> Evaluates the one dimensional source term of the density*velocity.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_1d_source_rho_u(x) bind (C,name='masa_eval_1d_source_rho_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       
     end function masa_eval_1d_source_rho_u
  end interface  

  interface 
     !> Evaluates the one dimensional source term of the energy*density.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_1d_source_rho_e(x) bind (C,name='masa_eval_1d_source_rho_e')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       
     end function masa_eval_1d_source_rho_e
  end interface

  interface 
     !> Evaluates the one dimensional source term of the 
     !! concentration of Nitrogen.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] funct Procedure that returns value of 
     !! equilibrium constant as a function of temperature, e.g. f(T)
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_1d_source_rho_N(x,funct) bind (C,name='masa_eval_1d_source_rho_N')
       use iso_c_binding
       implicit none
       
       real (c_double), value    :: x
       real (c_double), external :: funct
       
     end function masa_eval_1d_source_rho_N
  end interface  

  interface 
 
     !> Evaluates the one dimensional source term of the 
     !! concentration of Nitrogen-two.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] funct Procedure that returns value of 
     !! equilibrium constant as a function of temperature, e.g. f(T)
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_1d_source_rho_N2(x,funct) bind (C,name='masa_eval_1d_source_rho_N2')
       use iso_c_binding
       implicit none
       
       real (c_double), value    :: x
       real (c_double), external :: funct
       
     end function masa_eval_1d_source_rho_N2
  end interface  

  ! ---------------------------------
  ! MMS source term interfaces -- 2d
  ! ---------------------------------

  interface 
     !> Evaluates the two dimensional source term of the temperature.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_2d_source_t(x,y) bind (C,name='masa_eval_2d_source_t')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       
     end function masa_eval_2d_source_t
  end interface

  interface 
     !> Evaluates the two dimensional source term of the u-component
     !! of velocity.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_2d_source_u(x,y) bind (C,name='masa_eval_2d_source_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y

     end function masa_eval_2d_source_u
  end interface

  interface 
     !> Evaluates the two dimensional source term of the v-component
     !! of velocity.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_2d_source_v(x,y) bind (C,name='masa_eval_2d_source_v')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       
     end function masa_eval_2d_source_v
  end interface

  interface
     !> Evaluates the two dimensional source term of the energy.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_2d_source_e(x,y) bind (C,name='masa_eval_2d_source_e')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       
     end function masa_eval_2d_source_e
  end interface

  interface 
     !> Evaluates the two dimensional source term of the density.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_2d_source_rho(x,y) bind (C,name='masa_eval_2d_source_rho')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       
     end function masa_eval_2d_source_rho
  end interface  

  interface 
     !> Evaluates the two dimensional source term of the 
     !! density*u-component of velocity.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_2d_source_rho_u(x,y) bind (C,name='masa_eval_2d_source_rho_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       
     end function masa_eval_2d_source_rho_u
  end interface  

  interface 
     !> Evaluates the two dimensional source term of the 
     !! density*v-component of velocity.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_2d_source_rho_v(x,y) bind (C,name='masa_eval_2d_source_rho_v')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       
     end function masa_eval_2d_source_rho_v
  end interface
  
  interface 
     !> Evaluates the two dimensional source term of the 
     !! density*w-component of velocity.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_2d_source_rho_w(x,y) bind (C,name='masa_eval_2d_source_rho_w')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       
     end function masa_eval_2d_source_rho_w
  end interface  

  interface 
     !> Evaluates the two dimensional source term of the 
     !! density*energy.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_2d_source_rho_e(x,y) bind (C,name='masa_eval_2d_source_rho_e')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       
     end function masa_eval_2d_source_rho_e
  end interface  

  ! ---------------------------------
  ! MMS source term interfaces -- 3d
  ! ---------------------------------

  interface  
     !> Evaluates the three dimensional source term of the 
     !! temperature.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @param[in] z Real(8) value of the z-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_3d_source_t(x,y,z) bind (C,name='masa_eval_3d_source_t')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z
       
     end function masa_eval_3d_source_t
  end interface

  interface 
     !> Evaluates the three dimensional source term of the 
     !! u-component of velocity.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @param[in] z Real(8) value of the z-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_3d_source_u(x,y,z) bind (C,name='masa_eval_3d_source_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z
       
     end function masa_eval_3d_source_u
  end interface

  interface 
     !> Evaluates the three dimensional source term of the 
     !! v-component of velocity.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @param[in] z Real(8) value of the z-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_3d_source_v(x,y,z) bind (C,name='masa_eval_3d_source_v')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z

     end function masa_eval_3d_source_v
  end interface

  interface 
     !> Evaluates the three dimensional source term of the 
     !! w-component of velocity.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @param[in] z Real(8) value of the z-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_3d_source_w(x,y,z) bind (C,name='masa_eval_3d_source_w')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z

     end function masa_eval_3d_source_w
  end interface

  interface 
     !> Evaluates the three dimensional source term of the 
     !! energy.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @param[in] z Real(8) value of the z-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_3d_source_e(x,y,z) bind (C,name='masa_eval_3d_source_e')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z
       
     end function masa_eval_3d_source_e
  end interface

  interface 
     !> Evaluates the three dimensional source term of the 
     !! density.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @param[in] z Real(8) value of the z-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_3d_source_rho(x,y,z) bind (C,name='masa_eval_3d_source_rho')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z

     end function masa_eval_3d_source_rho
  end interface  

  interface 
     !> Evaluates the three dimensional source term of the 
     !! density*u-component of velocity.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @param[in] z Real(8) value of the z-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_3d_source_rho_u(x,y,z) bind (C,name='masa_eval_3d_source_rho_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z

     end function masa_eval_3d_source_rho_u
  end interface  

  interface 
     !> Evaluates the three dimensional source term of the 
     !! density*v-component of velocity.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @param[in] z Real(8) value of the z-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_3d_source_rho_v(x,y,z) bind (C,name='masa_eval_3d_source_rho_v')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z

     end function masa_eval_3d_source_rho_v
  end interface  

  interface 
     !> Evaluates the three dimensional source term of the 
     !! density*w-component of velocity.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @param[in] z Real(8) value of the z-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_3d_source_rho_w(x,y,z) bind (C,name='masa_eval_3d_source_rho_w')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z

     end function masa_eval_3d_source_rho_w
  end interface  

  interface 
     !> Evaluates the three dimensional source term of the 
     !! density*energy.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @param[in] z Real(8) value of the z-coordinate.
     !! @return Real(8) value for the source term.
     !!
     real (c_double) function masa_eval_3d_source_rho_e(x,y,z) bind (C,name='masa_eval_3d_source_rho_e')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z
       
     end function masa_eval_3d_source_rho_e
  end interface  

  ! ---------------------------------
  ! MMS analytical term interfaces -- 1d
  ! ---------------------------------

  interface 
     !> Evaluates the one dimensional exact solution of the 
     !! temperature.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @return Real(8) value for the exact solution.
     !!
     real (c_double) function masa_eval_1d_exact_t(x) bind (C,name='masa_eval_1d_exact_t')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       
     end function masa_eval_1d_exact_t
  end interface

  interface 
     !> Evaluates the one dimensional exact solution of the 
     !! u-component of velocity.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @return Real(8) value for the exact solution.
     !!
     real (c_double) function masa_eval_1d_exact_u(x) bind (C,name='masa_eval_1d_exact_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       
     end function masa_eval_1d_exact_u
  end interface

  interface 
     !> Evaluates the one dimensional exact solution of the 
     !! pressure.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @return Real(8) value for the exact solution.
     !!
     real (c_double) function masa_eval_1d_exact_p(x) bind (C,name='masa_eval_1d_exact_p')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       
     end function masa_eval_1d_exact_p
  end interface

  interface 
     !> Evaluates the one dimensional exact solution of the 
     !! density.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @return Real(8) value for the exact solution.
     !!
     real (c_double) function masa_eval_1d_exact_rho(x) bind (C,name='masa_eval_1d_exact_rho')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       
     end function masa_eval_1d_exact_rho
  end interface

  interface 
     !> Evaluates the one dimensional exact solution of the 
     !! Nitrogen.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @return Real(8) value for the exact solution.
     !!
     real (c_double) function masa_eval_1d_exact_rho_N(x) bind (C,name='masa_eval_1d_exact_rho_N')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       
     end function masa_eval_1d_exact_rho_N
  end interface

  interface 
     !> Evaluates the one dimensional exact solution of the 
     !! Nitrogen two.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @return Real(8) value for the exact solution.
     !!
     real (c_double) function masa_eval_1d_exact_rho_N2(x) bind (C,name='masa_eval_1d_exact_rho_N2')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       
     end function masa_eval_1d_exact_rho_N2
  end interface

  ! ---------------------------------
  ! MMS analytical term interfaces -- 2d
  ! ---------------------------------

  interface 
     !> Evaluates the two dimensional exact solution of the 
     !! temperature.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @return Real(8) value for the exact solution.
     !!
     real (c_double) function masa_eval_2d_exact_t(x,y) bind (C,name='masa_eval_2d_exact_t')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       
     end function masa_eval_2d_exact_t
  end interface

  interface 
     !> Evaluates the two dimensional exact solution of the 
     !! u-component of velocity.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @return Real(8) value for the exact solution.
     !!
     real (c_double) function masa_eval_2d_exact_u(x,y) bind (C,name='masa_eval_2d_exact_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       
     end function masa_eval_2d_exact_u
  end interface

  interface 
     !> Evaluates the two dimensional exact solution of the 
     !! v-component of velocity.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @return Real(8) value for the exact solution.
     !!
     real (c_double) function masa_eval_2d_exact_v(x,y) bind (C,name='masa_eval_2d_exact_v')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       
     end function masa_eval_2d_exact_v
  end interface

  interface 
     !> Evaluates the two dimensional exact solution of the 
     !! pressure.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @return Real(8) value for the exact solution.
     !!
     real (c_double) function masa_eval_2d_exact_p(x,y) bind (C,name='masa_eval_2d_exact_p')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       
     end function masa_eval_2d_exact_p
  end interface

  interface 
     !> Evaluates the two dimensional exact solution of the 
     !! density.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @return Real(8) value for the exact solution.
     !!
     real (c_double) function masa_eval_2d_exact_rho(x,y) bind (C,name='masa_eval_2d_exact_rho')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       
     end function masa_eval_2d_exact_rho
  end interface

  ! ---------------------------------
  ! MMS analytical term interfaces -- 3d
  ! ---------------------------------

  interface 
     !> Evaluates the three dimensional exact solution of the 
     !! temperature.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @param[in] z Real(8) value of the z-coordinate.
     !! @return Real(8) value for the exact solution.
     !!
     real (c_double) function masa_eval_3d_exact_t(x,y,z) bind (C,name='masa_eval_3d_exact_t')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z       

     end function masa_eval_3d_exact_t
  end interface

  interface 
     !> Evaluates the three dimensional exact solution of the 
     !! u-component of velocity.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @param[in] z Real(8) value of the z-coordinate.
     !! @return Real(8) value for the exact solution.
     !!
     real (c_double) function masa_eval_3d_exact_u(x,y,z) bind (C,name='masa_eval_3d_exact_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z       
       
     end function masa_eval_3d_exact_u
  end interface

  interface 
     !> Evaluates the three dimensional exact solution of the 
     !! v-component of velocity.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @param[in] z Real(8) value of the z-coordinate.
     !! @return Real(8) value for the exact solution.
     !!
     real (c_double) function masa_eval_3d_exact_v(x,y,z) bind (C,name='masa_eval_3d_exact_v')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z       
       
     end function masa_eval_3d_exact_v
  end interface

  interface 
     !> Evaluates the three dimensional exact solution of the 
     !! w-component of velocity.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @param[in] z Real(8) value of the z-coordinate.
     !! @return Real(8) value for the exact solution.
     !!
     real (c_double) function masa_eval_3d_exact_w(x,y,z) bind (C,name='masa_eval_3d_exact_w')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z       
       
     end function masa_eval_3d_exact_w
  end interface

  interface
     !> Evaluates the three dimensional exact solution of the 
     !! pressure.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @param[in] z Real(8) value of the z-coordinate.
     !! @return Real(8) value for the exact solution.
     !!
     real (c_double) function masa_eval_3d_exact_p(x,y,z) bind (C,name='masa_eval_3d_exact_p')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z       
       
     end function masa_eval_3d_exact_p
  end interface

  interface 
     !> Evaluates the three dimensional exact solution of the 
     !! density.
     !!
     !! @param[in] x Real(8) value of the x-coordinate.
     !! @param[in] y Real(8) value of the y-coordinate.
     !! @param[in] z Real(8) value of the z-coordinate.
     !! @return Real(8) value for the exact solution.
     !!
     real (c_double) function masa_eval_3d_exact_rho(x,y,z) bind (C,name='masa_eval_3d_exact_rho')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       real (c_double), value :: y
       real (c_double), value :: z       
       
     end function masa_eval_3d_exact_rho
  end interface

  ! ---------------------------------
  ! MMS gradient term interfaces -- 1d
  ! ---------------------------------

  interface 
     real (c_double) function masa_eval_1d_grad_u(x) bind (C,name='masa_eval_1d_grad_u')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       
     end function masa_eval_1d_grad_u
  end interface

  interface 
     real (c_double) function masa_eval_1d_grad_p(x) bind (C,name='masa_eval_1d_grad_p')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       
     end function masa_eval_1d_grad_p
  end interface

  interface 
     real (c_double) function masa_eval_1d_grad_rho(x) bind (C,name='masa_eval_1d_grad_rho')
       use iso_c_binding
       implicit none
       
       real (c_double), value :: x
       
     end function masa_eval_1d_grad_rho
  end interface

  ! ---------------------------------
  ! MMS gradient term interfaces -- 2d
  ! ---------------------------------

  interface 
     real (c_double) function masa_eval_2d_grad_u(x,y,it) bind (C,name='masa_eval_2d_grad_u')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: x
       real    (c_double), value :: y
       integer (c_int),    value :: it
       
     end function masa_eval_2d_grad_u
  end interface

  interface 
     real (c_double) function masa_eval_2d_grad_v(x,y,it) bind (C,name='masa_eval_2d_grad_v')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: x
       real    (c_double), value :: y
       integer (c_int),    value :: it
       
     end function masa_eval_2d_grad_v
  end interface

  interface 
     real (c_double) function masa_eval_2d_grad_w(x,y,it) bind (C,name='masa_eval_2d_grad_w')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: x
       real    (c_double), value :: y
       integer (c_int),    value :: it
       
     end function masa_eval_2d_grad_w
  end interface

  interface 
     real (c_double) function masa_eval_2d_grad_p(x,y,it) bind (C,name='masa_eval_2d_grad_p')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: x
       real    (c_double), value :: y
       integer (c_int),    value :: it
       
     end function masa_eval_2d_grad_p
  end interface

  interface 
     real (c_double) function masa_eval_2d_grad_rho(x,y,it) bind (C,name='masa_eval_2d_grad_rho')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: x
       real    (c_double), value :: y
       integer (c_int),    value :: it
       
     end function masa_eval_2d_grad_rho
  end interface

  ! ---------------------------------
  ! MMS gradient term interfaces -- 3d
  ! ---------------------------------

  interface 
     real (c_double) function masa_eval_3d_grad_u(x,y,z,it) bind (C,name='masa_eval_3d_grad_u')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: x
       real    (c_double), value :: y
       real    (c_double), value :: z
       integer (c_int),    value :: it
       
     end function masa_eval_3d_grad_u
  end interface

  interface 
     real (c_double) function masa_eval_3d_grad_v(x,y,z,it) bind (C,name='masa_eval_3d_grad_v')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: x
       real    (c_double), value :: y
       real    (c_double), value :: z
       integer (c_int),    value :: it
       
     end function masa_eval_3d_grad_v
  end interface

  interface 
     real (c_double) function masa_eval_3d_grad_w(x,y,z,it) bind (C,name='masa_eval_3d_grad_w')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: x
       real    (c_double), value :: y
       real    (c_double), value :: z
       integer (c_int),    value :: it
       
     end function masa_eval_3d_grad_w
  end interface

  interface 
     real (c_double) function masa_eval_3d_grad_p(x,y,z,it) bind (C,name='masa_eval_3d_grad_p')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: x
       real    (c_double), value :: y
       real    (c_double), value :: z
       integer (c_int),    value :: it
       
     end function masa_eval_3d_grad_p
  end interface

  interface 
     real (c_double) function masa_eval_3d_grad_rho(x,y,z,it) bind (C,name='masa_eval_3d_grad_rho')
       use iso_c_binding
       implicit none
       
       real    (c_double), value :: x
       real    (c_double), value :: y
       real    (c_double), value :: z
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

  !! \name sets the parameter value
  subroutine masa_set_param(param_name,value)
    use iso_c_binding
    implicit none

    character(len=*), intent(in)        :: param_name
    real  (c_double), intent(in)        :: value

    call masa_set_param_passthrough(param_name//C_NULL_CHAR,value)

  end subroutine masa_set_param

  ! ---------------------------------
  ! MMS Vector/Array Routines
  ! ---------------------------------

  subroutine masa_get_array(param_name,it,arr)
    use iso_c_binding
    implicit none
    
    character(len=*)                :: param_name
    integer (c_int)                 :: it          ! pass-by-ref is intentional
    real    (c_double),dimension(*) :: arr

    call masa_get_array_passthrough(param_name//C_NULL_CHAR, it, arr)
    
  end subroutine masa_get_array

end module masa
