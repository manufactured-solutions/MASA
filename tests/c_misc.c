// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
// $Author$
// $Id$
//
// c_misc.c :program that tests masa helper functions
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


int main()
{


  // list all currently initialized mms
  cmasa_list_mms();

  // initalize two different functions
  cmasa_init("euler-test","euler_1d");
  cmasa_init_param();
  cmasa_init("euler-test2","euler_2d");
  cmasa_init_param();

  cmasa_list_mms();

  // display parameters @ default values
  cmasa_display_param();
  
  //change parameter
  cmasa_set_param("u_0",2.3);

  // display again
  cmasa_display_param();

  //switch to original function
  cmasa_select_mms("euler-test");
  cmasa_display_param();

  //tests passed
  return 0;

}

