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
// This header file contains the public C-interface functions designed to be exposed in MASA
//

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
int cmasa_init_param();
//int cmasa_set_param(char*,double);
//int cmasa_get_param(char*,double*);

// --------------------------------
// source term(s) -- 1D
// --------------------------------
int cmasa_eval_t_source  (double,double*);        // x
int cmasa_eval_t_source  (double,double,double*); // x,t
int cmasa_eval_u_source  (double,double*);
int cmasa_eval_e_source  (double,double*);
int cmasa_eval_rho_source(double,double*);

int cmasa_eval_t_an      (double,double*);        // x
int cmasa_eval_t_an      (double,double,double*); // x,t
int cmasa_eval_u_an      (double,double*);
int cmasa_eval_p_an      (double,double*);
int cmasa_eval_rho_an    (double,double*);

// --------------------------------
// source term(s) -- 2D
// --------------------------------
int cmasa_eval_t_source  (double,double,double,double*); //x,y,t
int cmasa_eval_u_source  (double,double,double*);
int cmasa_eval_v_source  (double,double,double*);
int cmasa_eval_e_source  (double,double,double*);
int cmasa_eval_rho_source(double,double,double*);

int cmasa_eval_t_an      (double,double,double,double*); //x,y,t
int cmasa_eval_u_an      (double,double,double*);
int cmasa_eval_v_an      (double,double,double*);
int cmasa_eval_p_an      (double,double,double*);
int cmasa_eval_rho_an    (double,double,double*);

// --------------------------------
// source term(s) -- 3D
// --------------------------------
int cmasa_eval_t_source  (double,double,double,double,double*); // x,y,z,t
int cmasa_eval_u_source  (double,double,double,double*);
int cmasa_eval_v_source  (double,double,double,double*);
int cmasa_eval_w_source  (double,double,double,double*);
int cmasa_eval_e_source  (double,double,double,double*);
int cmasa_eval_rho_source(double,double,double,double*);

int cmasa_eval_t_an      (double,double,double,double,double*); // x,y,z,t
int cmasa_eval_u_an      (double,double,double,double*);
int cmasa_eval_v_an      (double,double,double,double*);
int cmasa_eval_w_an      (double,double,double,double*);
int cmasa_eval_p_an      (double,double,double,double*);
int cmasa_eval_rho_an    (double,double,double,double*);

// --------------------------------
// gradient of analytical solution
// --------------------------------
int cmasa_eval_1d_grad(int,double,double*);
int cmasa_eval_2d_grad(int,double,double,double*);
int cmasa_eval_3d_grad(int,double,double,double,double*);

// --------------------------------
// old masa core functions (to be called by user)
// --------------------------------

int cmasa_get_name(const char*);
int cmasa_get_dimension(int*);

int cmasa_display_param();
int cmasa_sanity_check();
