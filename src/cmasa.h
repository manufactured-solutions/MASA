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
  extern int cmasa_eval_1d_t_source  (double,double*);        // x
  extern int cmasa_eval_1d_u_source  (double,double*);
  extern int cmasa_eval_1d_e_source  (double,double*);
  extern int cmasa_eval_1d_rho_source(double,double*);

  extern int cmasa_eval_1d_t_an      (double,double*);        // x
  extern int cmasa_eval_1d_u_an      (double,double*);
  extern int cmasa_eval_1d_p_an      (double,double*);
  extern int cmasa_eval_1d_rho_an    (double,double*);

  //int cmasa_eval_t_source  (double,double,double*); // x,t
  //int cmasa_eval_t_an      (double,double,double*); // x,t

  // --------------------------------
  // source term(s) -- 2D
  // --------------------------------

  extern int cmasa_eval_2d_u_source  (double,double,double*);
  extern int cmasa_eval_2d_v_source  (double,double,double*);
  extern int cmasa_eval_2d_e_source  (double,double,double*);
  extern int cmasa_eval_2d_rho_source(double,double,double*);
  //int cmasa_eval__2dt_source  (double,double,double,double*); //x,y,t

  extern int cmasa_eval_2d_t_an      (double,double,double,double*); //x,y,t
  extern int cmasa_eval_2d_u_an      (double,double,double*);
  extern int cmasa_eval_2d_v_an      (double,double,double*);
  extern int cmasa_eval_2d_p_an      (double,double,double*);
  extern int cmasa_eval_2d_rho_an    (double,double,double*);

  // --------------------------------
  // source term(s) -- 3D
  // --------------------------------

  //int cmasa_eval_t_source  (double,double,double,double,double*); // x,y,z,t
  extern int cmasa_eval_3d_u_source  (double,double,double,double*);
  extern int cmasa_eval_3d_v_source  (double,double,double,double*);
  extern int cmasa_eval_3d_w_source  (double,double,double,double*);
  extern int cmasa_eval_3d_e_source  (double,double,double,double*);
  extern int cmasa_eval_3d_rho_source(double,double,double,double*);

  //int cmasa_eval_t_an      (double,double,double,double,double*); // x,y,z,t
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
