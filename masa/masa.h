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
// This header file contains the public functions designed to be exposed in MASA
//

#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <map>
#include <stdlib.h>

using namespace std;

namespace MASA
{
  // masa core functions (to be called by user)
  int masa_getid(void**,string);
  int masa_printid();
  int masa_set_param(void*,string,double);
  int masa_get_param(void*,string,double*);
  int masa_display_param(void*);
  int masa_get_name(void*,string*);
  int masa_get_dimension(void*,int*);
  int masa_sanity_check(void*);

  // source term(s) -- 1D
  int masa_eval_t_source  (void*,double,double*);
  int masa_eval_u_source  (void*,double,double*);
  int masa_eval_e_source  (void*,double,double*);
  int masa_eval_rho_source(void*,double,double*);

  // source term(s) -- 2D
  int masa_eval_t_source  (void*,double,double,double*);
  int masa_eval_u_source  (void*,double,double,double*);
  int masa_eval_v_source  (void*,double,double,double*);
  int masa_eval_w_source  (void*,double,double,double*);
  int masa_eval_e_source  (void*,double,double,double*);
  int masa_eval_rho_source(void*,double,double,double*);

  // source term(s) -- 3D
  int masa_eval_t_source  (void*,double,double,double,double*);
  int masa_eval_u_source  (void*,double,double,double,double*);
  int masa_eval_v_source  (void*,double,double,double,double*);
  int masa_eval_w_source  (void*,double,double,double,double*);
  int masa_eval_e_source  (void*,double,double,double,double*);
  int masa_eval_rho_source(void*,double,double,double,double*);

  // analytical solution
  int masa_eval_1d_an(void*,double,double*);
  int masa_eval_2d_an(void*,double,double,double*);
  int masa_eval_3d_an(void*,double,double,double,double*);

  // gradient of analytical solution
  int masa_eval_1d_grad(void*,int,double,double*);
  int masa_eval_2d_grad(void*,int,double,double,double*);
  int masa_eval_3d_grad(void*,int,double,double,double,double*);
  
} //end MASA namespace
