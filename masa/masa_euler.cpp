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
//   These are the MASA class member functions and constructors
//   For the Euler Equations

#include <masa_internal.h> 

using namespace MASA;

/* ------------------------------------------------
 *
 *         EULER EQUATION
 *
 *
 *
 * -----------------------------------------------
 */ 

MASA::euler_2d::euler_2d()
{
  mmsname = "euler_2d";
  dimension=2;

  register_var("u_0",&u_0);
  register_var("u_x",&u_x);
  register_var("u_y",&u_y);
  register_var("v_0",&v_0);
  register_var("v_x",&v_x);
  register_var("v_y",&v_y);
  register_var("rho_0",&rho_0);
  register_var("rho_x",&rho_x);
  register_var("rho_y",&rho_y);
  register_var("p_0",&p_0);
  register_var("p_x",&p_x);
  register_var("p_y",&p_y);
  register_var("a_px",&a_px);
  register_var("a_py",&a_py);
  register_var("a_rhox",&a_rhox);
  register_var("a_rhoy",&a_rhoy);
  register_var("a_ux",&a_ux);
  register_var("a_uy",&a_uy);
  register_var("a_vx",&a_vx);
  register_var("a_vy",&a_vy);
  register_var("L",&L);
  register_var("Gamma",&Gamma);
  register_var("mu",&mu);

}//done with constructor

MASA::euler_3d::euler_3d()
{
  mmsname = "euler_3d";
  dimension=3;

  register_var("u_0",&u_0);
  register_var("u_x",&u_x);
  register_var("u_y",&u_y);
  register_var("u_z",&u_z);
  register_var("v_0",&v_0);
  register_var("v_x",&v_x);
  register_var("v_y",&v_y);
  register_var("v_z",&v_z);
  register_var("rho_0",&rho_0);
  register_var("rho_x",&rho_x);
  register_var("rho_y",&rho_y);
  register_var("rho_z",&rho_z);
  register_var("p_0",&p_0);
  register_var("p_x",&p_x);
  register_var("p_y",&p_y);
  register_var("p_z",&p_z);
  register_var("a_px",&a_px);
  register_var("a_py",&a_py);
  register_var("a_rhox",&a_rhox);
  register_var("a_rhoy",&a_rhoy);
  register_var("a_ux",&a_ux);
  register_var("a_uy",&a_uy);
  register_var("a_uz",&a_uz);
  register_var("a_vx",&a_vx);
  register_var("a_vy",&a_vy);
  register_var("L",&L);
  register_var("Gamma",&Gamma);
  register_var("mu",&mu);
  register_var("a_pz",&a_pz);
  register_var("a_rhoz",&a_rhoz);
  register_var("a_vz",&a_vz);
  register_var("a_wz",&a_wz);
  register_var("a_wx",&a_wx);
  register_var("a_wy",&a_wy);
  register_var("w_x",&w_x);
  register_var("w_y",&w_y);
  register_var("w_z",&w_z);

}//done with constructor

/*double MASA::heateq_1d_steady_const::eval_q_u(double x)
{
  double qt = ax * ax * k0 * cos(ax * x);
  return qt;

}
*/
