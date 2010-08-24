// $License$
// $Author$
// $Id$ 

  /*--------------------------------------------------------------------------
  *--------------------------------------------------------------------------
  *
  * Copyright (C) 2010 The PECOS Development Team
  *
  * Please see http://pecos.ices.utexas.edu for more information.
  *
  * This file is part of MASA.
  *
  * MASA is free software: you can redistribute it and/or modify it under
  * the terms of the GNU Lesser General Public License as published by the Free
  * Software Foundation, either version 3 of the License, or (at your option)
  * any later version.
  *
  * MASA is distributed in the hope that it will be useful, but WITHOUT ANY
  * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
  * details.
  *
  * You should have received a copy of the GNU Lesser General Public License along 
  * with MASA.  If not, see <http://www.gnu.org/licenses/>.
  *
  *--------------------------------------------------------------------------
  
  MASA -- Manufactured Analytical Solutions Abstraction Library

  A software interface that provides access to all manufactured solutions to 
  be used by various models throughout the center.
  
  *--------------------------------------------------------------------------
  */  


/* 
  - this is an example of the heat equation in 1d for C-language -
*/


#include <cmasa.h>
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
  cmasa_get_param(a,&q);
  printf("A_x is set to: %g\n", q);

  // now lets change that parameters value to something else.
  cmasa_set_param(a,1.984);
  cmasa_get_param(a,&q);
  printf("A_x is set to: %g\n", q);

  //check all initialized properly
  cmasa_sanity_check();
  error = cmasa_eval_1d_t_source(x,&sol);
  printf("\nt source: %g\n",sol);
  
  cmasa_select_mms("bob");
  cmasa_display_param();
  error = cmasa_eval_1d_t_source(x2,&sol);
  printf("\nt source: %g\n",sol);
  return 0; // done
}
