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
// $Author: nick $
// $Id: masa_core.cpp 12639 2010-08-24 23:33:29Z nick $
//
// cmasa.h: contains the public C-interface functions in MASA
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#ifdef __cplusplus
extern "C" {
#endif
  
  // --------------------------------
  // new masa api function
  // --------------------------------
  extern int cmasa_init      (const char*, const char*);
  extern int cmasa_select_mms(const char*);
  extern int cmasa_curr_mms  (const char*);
  extern int cmasa_list_mms  ();

  // --------------------------------
  // interact with mms variables
  // --------------------------------
  extern int cmasa_init_param();
  extern int cmasa_set_param(const char*,double);
  extern int cmasa_get_param(const char*,double*);

  // --------------------------------
  // source term(s) -- 1D
  // --------------------------------
  extern int cmasa_eval_1d_t_source  (double,double*);
  extern int cmasa_eval_1d_u_source  (double,double*);
  extern int cmasa_eval_1d_e_source  (double,double*);
  extern int cmasa_eval_1d_rho_source(double,double*);

  extern int cmasa_eval_1d_t_an      (double,double*);
  extern int cmasa_eval_1d_u_an      (double,double*);
  extern int cmasa_eval_1d_p_an      (double,double*);
  extern int cmasa_eval_1d_rho_an    (double,double*);

  // --------------------------------
  // source term(s) -- 2D
  // --------------------------------

  extern int cmasa_eval_2d_u_source  (double,double,double*);
  extern int cmasa_eval_2d_v_source  (double,double,double*);
  extern int cmasa_eval_2d_e_source  (double,double,double*);
  extern int cmasa_eval_2d_rho_source(double,double,double*);

  extern int cmasa_eval_2d_t_an      (double,double,double*);
  extern int cmasa_eval_2d_u_an      (double,double,double*);
  extern int cmasa_eval_2d_v_an      (double,double,double*);
  extern int cmasa_eval_2d_p_an      (double,double,double*);
  extern int cmasa_eval_2d_rho_an    (double,double,double*);

  // --------------------------------
  // source term(s) -- 3D
  // --------------------------------

  extern int cmasa_eval_t_source     (double,double,double,double*);
  extern int cmasa_eval_3d_u_source  (double,double,double,double*);
  extern int cmasa_eval_3d_v_source  (double,double,double,double*);
  extern int cmasa_eval_3d_w_source  (double,double,double,double*);
  extern int cmasa_eval_3d_e_source  (double,double,double,double*);
  extern int cmasa_eval_3d_rho_source(double,double,double,double*);

  extern int cmasa_eval_t_an         (double,double,double,double*);
  extern int cmasa_eval_3d_u_an      (double,double,double,double*);
  extern int cmasa_eval_3d_v_an      (double,double,double,double*);
  extern int cmasa_eval_3d_w_an      (double,double,double,double*);
  extern int cmasa_eval_3d_p_an      (double,double,double,double*);
  extern int cmasa_eval_3d_rho_an    (double,double,double,double*);

  // --------------------------------
  // gradient of analytical solution
  // --------------------------------
  extern int cmasa_eval_1d_grad(int,double,double*);
  extern int cmasa_eval_2d_grad(int,double,double,double*);
  extern int cmasa_eval_3d_grad(int,double,double,double,double*);

  // --------------------------------
  // old masa core functions (to be called by user)
  // --------------------------------
  
  extern int cmasa_get_name(const char*);
  extern int cmasa_get_dimension(int*);
  extern int cmasa_display_param();
  extern int cmasa_sanity_check();

#ifdef __cplusplus
}
#endif
