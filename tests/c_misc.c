// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011 The PECOS Development Team
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

  // reroute stdout for regressions: TODO remove when logger mechanism
  // is used inside masa; these tests currently just verify functions
  // run successfully.

  freopen("/dev/null","w",stdout);

  // list all currently initialized mms
  masa_list_mms();

  // initalize two different functions
  masa_init("euler-test","euler_1d");
  masa_init_param();
  masa_init("euler-test2","euler_2d");
  masa_init_param();

  masa_list_mms();

  // display parameters @ default values
  masa_display_param();
  
  //change parameter
  masa_set_param("u_0",2.3);

  // display again
  masa_display_param();

  //switch to original function
  masa_select_mms("euler-test");
  masa_display_param();

  //tests passed
  return 0;

}

