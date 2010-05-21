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
//   For the Compressible Navier-Stokes

#include <masa_internal.h> 

using namespace MASA;

/* ------------------------------------------------
 *
 *         Compressible Navier Stokes Equations
 *
 *         2D
 *
 * -----------------------------------------------
 */ 

MASA::navierstokes_2d_compressible::navierstokes_2d_compressible()
{
  mmsname = "navierstokes_2d_compressible";
  dimension=2;

  //first variable (dummy) "dummy" -- load map and array
  dummy=MASA_VAR_DEFAULT;
  vararr.push_back(&dummy);

  // initalize other variables
  varmap["u_0"]=1;
  u_0=MASA_VAR_DEFAULT;              // need to initialize all variables!
  vararr.push_back(&u_0);

  // 2nd var
  varmap["u_x"]=2;
  u_x=MASA_VAR_DEFAULT;
  vararr.push_back(&u_x);

  // 3rd var
  varmap["u_y"]=3;
  u_y=MASA_VAR_DEFAULT;
  vararr.push_back(&u_y);

  // 4th var
  varmap["v_0"]=4;
  v_0=MASA_VAR_DEFAULT;
  vararr.push_back(&v_0);

  // 5th var
  varmap["v_x"]=5;
  v_x=MASA_VAR_DEFAULT;
  vararr.push_back(&v_x);

  //6th var
  varmap["v_y"]=6;
  v_y=MASA_VAR_DEFAULT;
  vararr.push_back(&v_y);

  //7th var
  varmap["rho_0"]=7;
  rho_0=MASA_VAR_DEFAULT;
  vararr.push_back(&rho_0);

  //8th var
  varmap["rho_x"]=8;
  rho_x=MASA_VAR_DEFAULT;
  vararr.push_back(&rho_x);

  //9th var
  varmap["rho_y"]=9;
  rho_y=MASA_VAR_DEFAULT;
  vararr.push_back(&rho_y);

  //10th var
  varmap["p_0"]=10;
  p_0=MASA_VAR_DEFAULT;
  vararr.push_back(&p_0);

  //11th var
  varmap["p_x"]=11;
  p_x=MASA_VAR_DEFAULT;
  vararr.push_back(&p_x);

  //12th var
  varmap["p_y"]=12;
  p_y=MASA_VAR_DEFAULT;
  vararr.push_back(&p_y);

  //13th var
  varmap["a_px"]=13;
  a_px=MASA_VAR_DEFAULT;
  vararr.push_back(&a_px);

  //14th var
  varmap["a_py"]=14;
  a_py=MASA_VAR_DEFAULT;
  vararr.push_back(&a_py);

  //15th var
  varmap["a_rhox"]=15;
  a_rhox=MASA_VAR_DEFAULT;
  vararr.push_back(&a_rhox);

  //16th var
  varmap["a_rhoy"]=16;
  a_rhoy=MASA_VAR_DEFAULT;
  vararr.push_back(&a_rhoy);

  //17th var
  varmap["a_ux"]=17;
  a_ux=MASA_VAR_DEFAULT;
  vararr.push_back(&a_ux);

  //18th var
  varmap["a_uy"]=18;
  a_uy=MASA_VAR_DEFAULT;
  vararr.push_back(&a_uy);

  //19th var
  varmap["a_vx"]=19;
  a_vx=MASA_VAR_DEFAULT;
  vararr.push_back(&a_vx);

  //20th var
  varmap["a_vy"]=20;
  a_vy=MASA_VAR_DEFAULT;
  vararr.push_back(&a_vy);

  //21th var
  varmap["L"]=21;
  L=MASA_VAR_DEFAULT;
  vararr.push_back(&L);

  //22nd var
  varmap["Gamma"]=22;
  Gamma=MASA_VAR_DEFAULT;
  vararr.push_back(&Gamma);

  //23rd var
  varmap["mu"]=23;
  mu=MASA_VAR_DEFAULT;
  vararr.push_back(&mu);

}//done with constructor


/* ------------------------------------------------
 *
 *         Compressible Navier Stokes Equations
 *
 *         3D
 *
 * -----------------------------------------------
 */ 

MASA::navierstokes_3d_compressible::navierstokes_3d_compressible()
{
  mmsname = "navierstokes_3d_compressible";
  dimension=3;

  //first variable (dummy) "dummy" -- load map and array
  dummy=MASA_VAR_DEFAULT;
  vararr.push_back(&dummy);

  // initalize other variables
  varmap["u_0"]=1;
  u_0=MASA_VAR_DEFAULT;              // need to initialize all variables!
  vararr.push_back(&u_0);

  // 2nd var
  varmap["u_x"]=2;
  u_x=MASA_VAR_DEFAULT;
  vararr.push_back(&u_x);

  // 3rd var
  varmap["u_y"]=3;
  u_y=MASA_VAR_DEFAULT;
  vararr.push_back(&u_y);

  // 4th var
  varmap["v_0"]=4;
  v_0=MASA_VAR_DEFAULT;
  vararr.push_back(&v_0);

  // 5th var
  varmap["v_x"]=5;
  v_x=MASA_VAR_DEFAULT;
  vararr.push_back(&v_x);

  //6th var
  varmap["v_y"]=6;
  v_y=MASA_VAR_DEFAULT;
  vararr.push_back(&v_y);

  //7th var
  varmap["rho_0"]=7;
  rho_0=MASA_VAR_DEFAULT;
  vararr.push_back(&rho_0);

  //8th var
  varmap["rho_x"]=8;
  rho_x=MASA_VAR_DEFAULT;
  vararr.push_back(&rho_x);

  //9th var
  varmap["rho_y"]=9;
  rho_y=MASA_VAR_DEFAULT;
  vararr.push_back(&rho_y);

  //10th var
  varmap["p_0"]=10;
  p_0=MASA_VAR_DEFAULT;
  vararr.push_back(&p_0);

  //11th var
  varmap["p_x"]=11;
  p_x=MASA_VAR_DEFAULT;
  vararr.push_back(&p_x);

  //12th var
  varmap["p_y"]=12;
  p_y=MASA_VAR_DEFAULT;
  vararr.push_back(&p_y);

  //13th var
  varmap["a_px"]=13;
  a_px=MASA_VAR_DEFAULT;
  vararr.push_back(&a_px);

  //14th var
  varmap["a_py"]=14;
  a_py=MASA_VAR_DEFAULT;
  vararr.push_back(&a_py);

  //15th var
  varmap["a_rhox"]=15;
  a_rhox=MASA_VAR_DEFAULT;
  vararr.push_back(&a_rhox);

  //16th var
  varmap["a_rhoy"]=16;
  a_rhoy=MASA_VAR_DEFAULT;
  vararr.push_back(&a_rhoy);

  //17th var
  varmap["a_ux"]=17;
  a_ux=MASA_VAR_DEFAULT;
  vararr.push_back(&a_ux);

  //18th var
  varmap["a_uy"]=18;
  a_uy=MASA_VAR_DEFAULT;
  vararr.push_back(&a_uy);

  //19th var
  varmap["a_vx"]=19;
  a_vx=MASA_VAR_DEFAULT;
  vararr.push_back(&a_vx);

  //20th var
  varmap["a_vy"]=20;
  a_vy=MASA_VAR_DEFAULT;
  vararr.push_back(&a_vy);

  //21th var
  varmap["L"]=21;
  L=MASA_VAR_DEFAULT;
  vararr.push_back(&L);

  //22th var
  varmap["u_z"]=22;
  u_z=MASA_VAR_DEFAULT;
  vararr.push_back(&u_z);

  //23th var
  varmap["v_z"]=23;
  v_z=MASA_VAR_DEFAULT;
  vararr.push_back(&v_z);

  //24th var
  varmap["w_z"]=24;
  w_z=MASA_VAR_DEFAULT;
  vararr.push_back(&w_z);

  //25th var
  varmap["rho_z"]=25;
  rho_z=MASA_VAR_DEFAULT;
  vararr.push_back(&rho_z);

  //26th var
  varmap["p_z"]=26;
  p_z=MASA_VAR_DEFAULT;
  vararr.push_back(&p_z);

  //27th var
  varmap["a_pz"]=27;
  a_pz=MASA_VAR_DEFAULT;
  vararr.push_back(&a_pz);

  //28th var
  varmap["a_rhoz"]=28;
  a_rhoz=MASA_VAR_DEFAULT;
  vararr.push_back(&a_rhoz);

  //29th var
  varmap["a_uz"]=29;
  a_uz=MASA_VAR_DEFAULT;
  vararr.push_back(&a_uz);

  //30th var
  varmap["a_vz"]=30;
  a_vz=MASA_VAR_DEFAULT;
  vararr.push_back(&a_vz);

  //31th var
  varmap["a_wz"]=31;
  a_wz=MASA_VAR_DEFAULT;
  vararr.push_back(&a_wz);

  //32nd var
  varmap["Gamma"]=32;
  Gamma=MASA_VAR_DEFAULT;
  vararr.push_back(&Gamma);

  //33rd var
  varmap["mu"]=33;
  mu=MASA_VAR_DEFAULT;
  vararr.push_back(&mu);

  //34th var
  varmap["a_wy"]=34;
  a_wy=MASA_VAR_DEFAULT;
  vararr.push_back(&a_wy);

  //35th var
  varmap["a_wx"]=35;
  a_wx=MASA_VAR_DEFAULT;
  vararr.push_back(&a_wx);

  //36th var
  varmap["a_py"]=36;
  a_py=MASA_VAR_DEFAULT;
  vararr.push_back(&a_py);

  //37th var
  varmap["a_px"]=37;
  a_px=MASA_VAR_DEFAULT;
  vararr.push_back(&a_px);

}//done with constructor

/*double MASA::heateq_1d_steady_const::eval_q_u(double x)
{
  double qt = ax * ax * k0 * cos(ax * x);
  return qt;

}
*/
