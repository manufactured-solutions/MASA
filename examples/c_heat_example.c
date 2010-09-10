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
// c_heat_example.c: this is an example of the heat equation in 1d for C
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa.h>
#include <stdio.h>

int main()
{
  double sol=0;
  double q=0;
  int error;
  double x  = 1.2;
  double x2 = 2.2;
  char* a= "A_x";
  char* b= "k_0";
  
  sol = 0;

  // init
  cmasa_init("nick","heateq_1d_steady_const");
  cmasa_init_param();

  cmasa_init("bob","heateq_1d_steady_const");
  cmasa_init_param();
  
  // list
  cmasa_list_mms();

  // switch
  cmasa_select_mms("nick");
  cmasa_display_param();

  // lets examine a particular parameter 
  q=cmasa_get_param(a);
  printf("A_x is set to: %g\n", q);

  // now lets change that parameters value to something else.
  cmasa_set_param(a,1.984);
  q=cmasa_get_param(a);
  printf("A_x is set to: %g\n", q);

  //check all initialized properly
  cmasa_sanity_check();
  sol = cmasa_eval_1d_t_source(x);
  printf("\nt source: %g\n",sol);
  
  cmasa_select_mms("bob");
  cmasa_display_param();
  sol = cmasa_eval_1d_t_source(x2);
  printf("\nt source: %g\n",sol);
  return 0; // done
}
