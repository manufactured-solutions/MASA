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
// 

#include <masa_internal.h>

using namespace MASA;

void MASA::manufactured_solution::get_var(string var, double* sol)
{
  int selector;
  selector = varmap[var];    // find location in pointer array

    if(selector==0)
    {
      cout << "\nMASA ERROR: No such variable exists\n";
      // exit(1); this isnt really a fatal error
    }
  else 
    {
      *sol = *vararr[selector];   // set to value 
    } 

}// done with get_var function

void MASA::manufactured_solution::set_var(string var, double val)
{
  int selector;
  selector = varmap[var];    // find location in pointer array

  if(selector==0)
    {
      cout << "\nMASA ERROR: No such variable to be set\n";
      // exit(1);not really a fatal error
    }
  else 
    {
      *vararr[selector] = val;   // set variable to new value    
    } 
}// done with set_var function

//test problem -- not a real class!
MASA::MASA_Test::MASA_Test()
{
  // here, we load up the map so we can key to specific variables
  // using input
  mmsname = "MASA example function";
  
  //first variable "axp" -- load map and array
  varmap["axp"]=1;
  axp=-1;
  vararr.push_back(&axp);
  vararr.push_back(&axp);
  
}//done with constructor

/* ------------------------------------------------
 *
 *         HEAT EQUATION
 *
 *
 *
 * -----------------------------------------------
 */ 

MASA::heateq_1d_steady_const::heateq_1d_steady_const()
{
    mmsname = "heateq_1d_steady_const";

  //first variable "axp" -- load map and array
  varmap["axp"]=1;
  axp=-1;
  vararr.push_back(&axp);
  vararr.push_back(&axp);

  // initalize other variables
  varmap["ax"]=2;
  ax=-1;
  vararr.push_back(&ax);

  // 2nd var
  varmap["k0"]=3;
  k0=-1;
  vararr.push_back(&k0);   

}//done with constructor

double MASA::heateq_1d_steady_const::eval_q_u(double x)
{
  double qt = ax * ax * k0 * cos(ax * x);
  return qt;

}//done with constructor


MASA::heateq_2d_steady_const::heateq_2d_steady_const()
{
    mmsname = "heateq_2d_steady_const";


}//done with constructor


MASA::heateq_3d_steady_const::heateq_3d_steady_const()
{
    mmsname = "heateq_3d_steady_const";


}//done with constructor


/* ------------------------------------------------
 *
 *         HEAT EQUATION
 *
 *   Unsteady         
 *
 * -----------------------------------------------
 */ 

MASA::heateq_1d_unsteady_const::heateq_1d_unsteady_const()
{
    mmsname = "heateq_1d_unsteady_const";


}//done with constructor


MASA::heateq_2d_unsteady_const::heateq_2d_unsteady_const()
{
    mmsname = "heateq_2d_unsteady_const";


}//done with constructor


MASA::heateq_3d_unsteady_const::heateq_3d_unsteady_const()
{
    mmsname = "heateq_3d_unsteady_const";


}//done with constructor


MASA::heateq_1d_unsteady_var::heateq_1d_unsteady_var()
{
    mmsname = "heateq_1d_unsteady_var";


}//done with constructor


MASA::heateq_2d_unsteady_var::heateq_2d_unsteady_var()
{
    mmsname = "heateq_2d_unsteady_var";


}//done with constructor


MASA::heateq_3d_unsteady_var::heateq_3d_unsteady_var()
{
    mmsname = "heateq_3d_unsteady_var";


}//done with constructor


/* ------------------------------------------------
 *
 *         HEAT EQUATION
 *
 *        steady
 *
 *        variable coefficients
 *
 * -----------------------------------------------
 */ 


MASA::heateq_1d_steady_var::heateq_1d_steady_var()
{
    mmsname = "heateq_1d_steady_var";


}//done with constructor


MASA::heateq_2d_steady_var::heateq_2d_steady_var()
{
    mmsname = "heateq_2d_steady_var";


}//done with constructor


MASA::heateq_3d_steady_var::heateq_3d_steady_var()
{
    mmsname = "heateq_3d_steady_var";


}//done with constructor

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


}//done with constructor

MASA::euler_3d::euler_3d()
{
    mmsname = "euler_3d";


}//done with constructor

/* ------------------------------------------------
 *
 *         Compressible Navier-Stokes EQUATIONs
 *
 *
 *
 * -----------------------------------------------
 */ 

MASA::ns_compress_2d::ns_compress_2d()
{
  mmsname = "navierstokes_compressible_2d";


}//done with constructor

MASA::ns_compress_3d::ns_compress_3d()
{
  mmsname = "navierstokes_compressible_3d";


}//done with constructor
