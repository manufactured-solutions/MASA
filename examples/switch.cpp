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


// 
//   This C++ example demonstrates how to switch between two different solutions
//

// header contains all subroutine information
#include <masa.h>

// MASA was designed with a custom namespace
using namespace MASA;



using namespace std;

int main()
{
  double solution;

  // initialize first solution
  masa_init("alice","heateq_1d_steady_const");

  // initialize 2nd solution
  masa_init("bob"  ,"euler_2d");
  
  masa_select_mms("alice");
  masa_init_param();
  masa_display_param();
  masa_eval_t_source(1.2,&solution);
  cout << solution << endl;

  masa_select_mms("bob");
  masa_init_param();
  masa_display_param();
  masa_eval_u_source(1,1,&solution);
  cout << solution << endl;

}

