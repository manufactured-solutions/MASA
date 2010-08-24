
// $License$

// $Id$


//-*-c++-*-
//
//----------------------------------------------------------------begin-lic-
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
//--------------------------------------------------------------------------
//
// masa.h: public functions designed to be exposed in MASA
//
// $Id$
//--------------------------------------------------------------------------
//------------------------------------------------------------------end-lic-
  

// This header file contains the public functions designed to be exposed in MASA
// What follows is the masa.h doxygen documentation headers

/*! \file masa.h
\brief MASA header file containing all public C++ API

MASA.h is a header file that contains all the public objects and member functions for the c++ interfaces.

*/

/*! \fn int masa_init(string unique_name, string mms)
\brief Initialize a MASA manufactured solution class
\param unique_name This is a string that provides a unique identity for the initialized manufactured class, e.g. "bob"
\param mms This is the manufactured class to be initialized, e.g. "euler_1d"
*/

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

  // --------------------------------
  // new masa api function
  // --------------------------------
  int masa_init      (string, string);
  int masa_select_mms(string);
  int masa_curr_mms  (string*);
  int masa_list_mms  ();

  // --------------------------------
  // interact with mms variables
  // --------------------------------
  int masa_init_param();
  int masa_set_param(string,double);
  int masa_get_param(string,double*);

  // --------------------------------
  // source term(s) -- 1D
  // --------------------------------
  int masa_eval_t_source  (double,double*);        // x
  int masa_eval_t_source  (double,double,double*); // x,t
  int masa_eval_u_source  (double,double*);
  int masa_eval_e_source  (double,double*);
  int masa_eval_rho_source(double,double*);
  int masa_eval_rho_u_source(double,double,double*); // 1d sod: x,t

  int masa_eval_t_an      (double,double*);        // x
  int masa_eval_t_an      (double,double,double*); // x,t
  int masa_eval_u_an      (double,double*);
  int masa_eval_p_an      (double,double*);
  int masa_eval_rho_an    (double,double*);

  // --------------------------------
  // source term(s) -- 2D
  // --------------------------------
  int masa_eval_t_source    (double,double,double,double*); //x,y,t
  int masa_eval_u_source    (double,double,double*);
  int masa_eval_v_source    (double,double,double*);
  int masa_eval_e_source    (double,double,double*);
  int masa_eval_rho_source  (double,double,double*);

  int masa_eval_t_an        (double,double,double,double*); //x,y,t
  int masa_eval_u_an        (double,double,double*);
  int masa_eval_v_an        (double,double,double*);
  int masa_eval_p_an        (double,double,double*);
  int masa_eval_rho_an      (double,double,double*);

  // --------------------------------
  // source term(s) -- 3D
  // --------------------------------
  int masa_eval_t_source  (double,double,double,double,double*); // x,y,z,t
  int masa_eval_u_source  (double,double,double,double*);
  int masa_eval_v_source  (double,double,double,double*);
  int masa_eval_w_source  (double,double,double,double*);
  int masa_eval_e_source  (double,double,double,double*);
  int masa_eval_rho_source(double,double,double,double*);

  int masa_eval_t_an      (double,double,double,double,double*); // x,y,z,t
  int masa_eval_u_an      (double,double,double,double*);
  int masa_eval_v_an      (double,double,double,double*);
  int masa_eval_w_an      (double,double,double,double*);
  int masa_eval_p_an      (double,double,double,double*);
  int masa_eval_rho_an    (double,double,double,double*);

  // --------------------------------
  // gradient of analytical solution
  // --------------------------------
  int masa_eval_1d_grad(int,double,double*);
  int masa_eval_2d_grad(int,double,double,double*);
  int masa_eval_3d_grad(int,double,double,double,double*);

  // --------------------------------
  // internal masa functions user might want to call
  // --------------------------------
  int masa_map (string*);
  int masa_map2(string, string);

  // --------------------------------
  // old masa core functions (to be called by user)
  // --------------------------------
  int masa_getid(void**,string);
  int masa_printid();
  int masa_display_param();
  int masa_get_name(string*);
  int masa_get_dimension(int*);
  int masa_sanity_check();
  
} //end MASA namespace
