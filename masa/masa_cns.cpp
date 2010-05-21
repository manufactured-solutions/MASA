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

MASA::compressible_navierstokes_2d::compressible_navierstokes_2d()
{
  mmsname = "compressible_navierstokes_2d";
  dimension=2;

  //first variable (dummy) "axp" -- load map and array
  axp=MASA_VAR_DEFAULT;
  vararr.push_back(&axp);

  // initalize other variables
  varmap["u0"]=1;
  u0=MASA_VAR_DEFAULT;              // need to initialize all variables!
  vararr.push_back(&u0);

  // 2nd var
  varmap["ux"]=2;
  ux=MASA_VAR_DEFAULT;
  vararr.push_back(&ux);

  // 3rd var
  varmap["uy"]=3;
  uy=MASA_VAR_DEFAULT;
  vararr.push_back(&uy);

  // 4th var
  varmap["v0"]=4;
  v0=MASA_VAR_DEFAULT;
  vararr.push_back(&v0);

  // 5th var
  varmap["vx"]=5;
  vx=MASA_VAR_DEFAULT;
  vararr.push_back(&vx);

  //6th var
  varmap["vy"]=6;
  vy=MASA_VAR_DEFAULT;
  vararr.push_back(&vy);

  //7th var
  varmap["rho0"]=7;
  rho0=MASA_VAR_DEFAULT;
  vararr.push_back(&rho0);

  //8th var
  varmap["rhox"]=8;
  rhox=MASA_VAR_DEFAULT;
  vararr.push_back(&rhox);

  //9th var
  varmap["rhoy"]=9;
  rhoy=MASA_VAR_DEFAULT;
  vararr.push_back(&rhoy);

  //10th var
  varmap["p0"]=10;
  p0=MASA_VAR_DEFAULT;
  vararr.push_back(&p0);

  //11th var
  varmap["px"]=11;
  px=MASA_VAR_DEFAULT;
  vararr.push_back(&px);

  //12th var
  varmap["py"]=12;
  py=MASA_VAR_DEFAULT;
  vararr.push_back(&py);

  //13th var
  varmap["apx"]=13;
  apx=MASA_VAR_DEFAULT;
  vararr.push_back(&apx);

  //14th var
  varmap["apy"]=14;
  apy=MASA_VAR_DEFAULT;
  vararr.push_back(&apy);

  //15th var
  varmap["arhox"]=15;
  arhox=MASA_VAR_DEFAULT;
  vararr.push_back(&arhox);

  //16th var
  varmap["arhoy"]=16;
  arhoy=MASA_VAR_DEFAULT;
  vararr.push_back(&arhoy);

  //17th var
  varmap["aux"]=17;
  aux=MASA_VAR_DEFAULT;
  vararr.push_back(&aux);

  //18th var
  varmap["auy"]=18;
  auy=MASA_VAR_DEFAULT;
  vararr.push_back(&auy);

  //19th var
  varmap["avx"]=19;
  avx=MASA_VAR_DEFAULT;
  vararr.push_back(&avx);

  //20th var
  varmap["avy"]=20;
  avy=MASA_VAR_DEFAULT;
  vararr.push_back(&avy);

  //21st var
  varmap["l"]=21;
  l=MASA_VAR_DEFAULT;
  vararr.push_back(&l);

  //22nd var
  varmap["gamma"]=22;
  gamma=MASA_VAR_DEFAULT;
  vararr.push_back(&gamma);

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

MASA::compressible_navierstokes_3d::compressible_navierstokes_3d()
{
  mmsname = "compressible_navierstokes_3d";
  dimension=3;

  //first variable (dummy) "axp" -- load map and array
  axp=MASA_VAR_DEFAULT;
  vararr.push_back(&axp);

  // initalize other variables
  varmap["u0"]=1;
  u0=MASA_VAR_DEFAULT;              // need to initialize all variables!
  vararr.push_back(&u0);

  // 2nd var
  varmap["ux"]=2;
  ux=MASA_VAR_DEFAULT;
  vararr.push_back(&ux);

  // 3rd var
  varmap["uy"]=3;
  uy=MASA_VAR_DEFAULT;
  vararr.push_back(&uy);

  // 4th var
  varmap["v0"]=4;
  v0=MASA_VAR_DEFAULT;
  vararr.push_back(&v0);

  // 5th var
  varmap["vx"]=5;
  vx=MASA_VAR_DEFAULT;
  vararr.push_back(&vx);

  //6th var
  varmap["vy"]=6;
  vy=MASA_VAR_DEFAULT;
  vararr.push_back(&vy);

  //7th var
  varmap["rho0"]=7;
  rho0=MASA_VAR_DEFAULT;
  vararr.push_back(&rho0);

  //8th var
  varmap["rhox"]=8;
  rhox=MASA_VAR_DEFAULT;
  vararr.push_back(&rhox);

  //9th var
  varmap["rhoy"]=9;
  rhoy=MASA_VAR_DEFAULT;
  vararr.push_back(&rhoy);

  //10th var
  varmap["p0"]=10;
  p0=MASA_VAR_DEFAULT;
  vararr.push_back(&p0);

  //11th var
  varmap["px"]=11;
  px=MASA_VAR_DEFAULT;
  vararr.push_back(&px);

  //12th var
  varmap["py"]=12;
  py=MASA_VAR_DEFAULT;
  vararr.push_back(&py);

  //13th var
  varmap["apx"]=13;
  apx=MASA_VAR_DEFAULT;
  vararr.push_back(&apx);

  //14th var
  varmap["apy"]=14;
  apy=MASA_VAR_DEFAULT;
  vararr.push_back(&apy);

  //15th var
  varmap["arhox"]=15;
  arhox=MASA_VAR_DEFAULT;
  vararr.push_back(&arhox);

  //16th var
  varmap["arhoy"]=16;
  arhoy=MASA_VAR_DEFAULT;
  vararr.push_back(&arhoy);

  //17th var
  varmap["aux"]=17;
  aux=MASA_VAR_DEFAULT;
  vararr.push_back(&aux);

  //18th var
  varmap["auy"]=18;
  auy=MASA_VAR_DEFAULT;
  vararr.push_back(&auy);

  //19th var
  varmap["avx"]=19;
  avx=MASA_VAR_DEFAULT;
  vararr.push_back(&avx);

  //20th var
  varmap["avy"]=20;
  avy=MASA_VAR_DEFAULT;
  vararr.push_back(&avy);

  //21th var
  varmap["l"]=21;
  l=MASA_VAR_DEFAULT;
  vararr.push_back(&l);

  //22th var
  varmap["uz"]=22;
  uz=MASA_VAR_DEFAULT;
  vararr.push_back(&uz);

  //23th var
  varmap["vz"]=23;
  vz=MASA_VAR_DEFAULT;
  vararr.push_back(&vz);

  //24th var
  varmap["wz"]=24;
  wz=MASA_VAR_DEFAULT;
  vararr.push_back(&wz);

  //25th var
  varmap["rhoz"]=25;
  rhoz=MASA_VAR_DEFAULT;
  vararr.push_back(&rhoz);

  //26th var
  varmap["pz"]=26;
  pz=MASA_VAR_DEFAULT;
  vararr.push_back(&pz);

  //27th var
  varmap["apz"]=27;
  apz=MASA_VAR_DEFAULT;
  vararr.push_back(&apz);

  //28th var
  varmap["arhoz"]=28;
  arhoz=MASA_VAR_DEFAULT;
  vararr.push_back(&arhoz);

  //29th var
  varmap["auz"]=29;
  auz=MASA_VAR_DEFAULT;
  vararr.push_back(&auz);

  //30th var
  varmap["avz"]=30;
  avz=MASA_VAR_DEFAULT;
  vararr.push_back(&avz);

  //31th var
  varmap["awz"]=31;
  awz=MASA_VAR_DEFAULT;
  vararr.push_back(&awz);

  //32nd var
  varmap["gamma"]=32;
  gamma=MASA_VAR_DEFAULT;
  vararr.push_back(&gamma);

  //33rd var
  varmap["mu"]=33;
  mu=MASA_VAR_DEFAULT;
  vararr.push_back(&mu);

}//done with constructor

/*double MASA::heateq_1d_steady_const::eval_q_u(double x)
{
  double qt = ax * ax * k0 * cos(ax * x);
  return qt;

}
*/
