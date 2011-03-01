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
#include <stdlib.h>


double MASA_VAR_DEFAULT = -12345.67;
double uninit = -1.33;

void test(double input)
{
  if(input == MASA_VAR_DEFAULT)
    {
      exit(1);
    }

  if(input == uninit)
    {
      exit(1);
    }

}

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
  masa_init_param();

  masa_init("bob","heateq_1d_steady_const");
  masa_init_param();
  
  // list
  masa_list_mms();

  // switch
  masa_select_mms("nick");
  masa_display_param();

  // lets examine a particular parameter 
  q=masa_get_param(a);
  printf("A_x is set to: %g\n", q);

  // now lets change that parameters value to something else.
  masa_set_param(a,1.984);
  q=masa_get_param(a);
  printf("A_x is set to: %g\n", q);

  //check all initialized properly
  masa_sanity_check();
  sol = masa_eval_1d_source_t(x);
  printf("\nt source: %g\n",sol);
  
  masa_select_mms("bob");
  masa_display_param();
  sol = masa_eval_1d_source_t(x2);
  printf("\nt source: %g\n",sol);
  return 0; // done
}
