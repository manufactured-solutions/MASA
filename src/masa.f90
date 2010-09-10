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
     subroutine masa_select_mms(desired_mms_function) bind (C,name='cmasa_select_mms')
       use iso_c_binding
       implicit none

       character(c_char), dimension(*), intent(in) :: desired_mms_function

     end subroutine masa_select_mms
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

     subroutine masa_set_param_passthrough(param_name,value) bind (C,name='cmasa_set_param')
       use iso_c_binding
       implicit none

       character(c_char), intent(in) :: param_name(*)
       real (c_double), value        :: value

     end subroutine masa_set_param_passthrough
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
