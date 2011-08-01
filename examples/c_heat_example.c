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
// c_heat_example.c: this is an example of the heat equation in 1d for C
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa.h>

typedef double Scalar;

int main()
{
  double sol=0;
  double q=0;
  double x  = 1.2;
  double x2 = 2.2;
  char* a= "A_x";
  
  sol = 0;

  // init
  masa_init("nick","heateq_1d_steady_const");
  masa_init("bob","heateq_1d_steady_const");

  // reroute stdout: comment to display to screen
  freopen("/dev/null","w",stdout);
  
  // we can list initialized mms with 
   masa_list_mms();

  // switch
  masa_select_mms("nick");
  
  // we can display the parameter list with
  masa_display_param();

  // lets examine a particular parameter 
  q=masa_get_param(a);

  // now lets change that parameters value to something else.
  masa_set_param(a,1.984);
  q=masa_get_param(a);

  //check all initialized properly
  masa_sanity_check();
  sol = masa_eval_1d_source_t(x);
  
  masa_select_mms("bob");
  sol = masa_eval_1d_source_t(x2);
  return 0; // done
}
