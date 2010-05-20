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
//   For the HEAT EQUATION

#include <masa_internal.h> 

using namespace MASA;

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
    dimension=1;

  // this bears some explaination -- to make the index match between the 
  // map and the pointer array, the index must both be starting at _1_,
  // not zero, as is typical for c. Thus, we have to add a dummy variable here

  //first variable (dummy) "axp" -- load map and array
  axp=MASA_VAR_DEFAULT;
  vararr.push_back(&axp);

  // initalize other variables
  varmap["ax"]=1;
  ax=MASA_VAR_DEFAULT;              // need to initialize all variables!
  vararr.push_back(&ax);

  // 2nd var
  varmap["k0"]=2;
  k0=MASA_VAR_DEFAULT;
  vararr.push_back(&k0);   

}//done with constructor

double MASA::heateq_1d_steady_const::eval_q_u(double x)
{
  double qt = ax * ax * k0 * cos(ax * x);
  return qt;

}

MASA::heateq_2d_steady_const::heateq_2d_steady_const()
{
    mmsname = "heateq_2d_steady_const";
    dimension=2;

    //first variable (dummy) "axp" -- load map and array
    axp=MASA_VAR_DEFAULT;
    vararr.push_back(&axp);
    
    // initalize other variables
    varmap["ax"]=1;                   // incriment map location
    ax=MASA_VAR_DEFAULT;              // need to initialize all variables!
    vararr.push_back(&ax);            // add variable to pointer array
    
    // 2nd var
    varmap["k0"]=2;
    k0=MASA_VAR_DEFAULT;
    vararr.push_back(&k0);   
    
    // 3rd var
    varmap["by"]=3;
    by=MASA_VAR_DEFAULT;
    vararr.push_back(&by);   

}//done with constructor

double MASA::heateq_2d_steady_const::eval_q_u(double x)
{
 
  //Q_T = k_0 * cos(A_x * x) * cos(B_y * y) * (A_x * A_x + B_y * B_y);
  //T_an = cos(A_x * x) * cos(B_y * y);
  //gradT_an[0] = -A_x * cos(B_y * y) * sin(A_x * x);
  //gradT_an[1] = -B_y * cos(A_x * x) * sin(B_y * y);

  //double qt = ax * ax * k0 * cos(ax * x);
  return 1;

}


MASA::heateq_3d_steady_const::heateq_3d_steady_const()
{
    mmsname = "heateq_3d_steady_const";
    dimension=3;

    //first variable (dummy) "axp" -- load map and array
    axp=MASA_VAR_DEFAULT;
    vararr.push_back(&axp);
    
    // initalize other variables
    varmap["ax"]=1;                   // incriment map location
    ax=MASA_VAR_DEFAULT;              // need to initialize all variables!
    vararr.push_back(&ax);            // add variable to pointer array
    
    // 2nd var
    varmap["k0"]=2;
    k0=MASA_VAR_DEFAULT;
    vararr.push_back(&k0);   
    
    // 3rd var
    varmap["by"]=3;
    by=MASA_VAR_DEFAULT;
    vararr.push_back(&by);   

    // 4th var
    varmap["cz"]=4;
    cz=MASA_VAR_DEFAULT;
    vararr.push_back(&cz);   

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
  dimension=1;

  //first variable (dummy) "axp" -- load map and array
  axp=MASA_VAR_DEFAULT;
  vararr.push_back(&axp);
  
  // initalize other variables
  varmap["ax"]=1;                   // incriment map location
  ax=MASA_VAR_DEFAULT;              // need to initialize all variables!
  vararr.push_back(&ax);            // add variable to pointer array
  
  // 2nd var
  varmap["k0"]=2;
  k0=MASA_VAR_DEFAULT;
  vararr.push_back(&k0);   
  
  // 3rd var
  varmap["dt"]=3;
  dt=MASA_VAR_DEFAULT;
  vararr.push_back(&dt);   
  
  // 4th var
  varmap["cp0"]=4;
  cp0=MASA_VAR_DEFAULT;
  vararr.push_back(&cp0);   

  // 5th var
  varmap["at"]=5;
  at=MASA_VAR_DEFAULT;
  vararr.push_back(&at);

  // 6th var
  varmap["rho"]=6;
  rho=MASA_VAR_DEFAULT;
  vararr.push_back(&rho);

}//done with constructor


MASA::heateq_2d_unsteady_const::heateq_2d_unsteady_const()
{
  mmsname = "heateq_2d_unsteady_const";
  dimension=2;
    
  //first variable (dummy) "axp" -- load map and array
  axp=MASA_VAR_DEFAULT;
  vararr.push_back(&axp);

  // initalize other variables
  varmap["ax"]=1;                   // incriment map location
  ax=MASA_VAR_DEFAULT;              // need to initialize all variables!
  vararr.push_back(&ax);            // add variable to pointer array
  
  // 2nd var
  varmap["k0"]=2;
  k0=MASA_VAR_DEFAULT;
  vararr.push_back(&k0);   
  
  // 3rd var
  varmap["dt"]=3;
  dt=MASA_VAR_DEFAULT;
  vararr.push_back(&dt);   
  
  // 4th var
  varmap["cp0"]=4;
  cp0=MASA_VAR_DEFAULT;
  vararr.push_back(&cp0);   

  // 5th var
  varmap["at"]=5;
  at=MASA_VAR_DEFAULT;
  vararr.push_back(&at);

  // 6th var
  varmap["rho"]=6;
  rho=MASA_VAR_DEFAULT;
  vararr.push_back(&rho);

  // 7th var
  varmap["by"]=7;
  by=MASA_VAR_DEFAULT;
  vararr.push_back(&by);

  // 8th var
  varmap["bt"]=7;
  bt=MASA_VAR_DEFAULT;
  vararr.push_back(&bt);

}//done with constructor


MASA::heateq_3d_unsteady_const::heateq_3d_unsteady_const()
{
    mmsname = "heateq_3d_unsteady_const";
    dimension=3;
    
  //first variable (dummy) "axp" -- load map and array
  axp=MASA_VAR_DEFAULT;
  vararr.push_back(&axp);

  // initalize other variables
  varmap["ax"]=1;                   // incriment map location
  ax=MASA_VAR_DEFAULT;              // need to initialize all variables!
  vararr.push_back(&ax);            // add variable to pointer array
  
  // 2nd var
  varmap["k0"]=2;
  k0=MASA_VAR_DEFAULT;
  vararr.push_back(&k0);   
  
  // 3rd var
  varmap["dt"]=3;
  dt=MASA_VAR_DEFAULT;
  vararr.push_back(&dt);   
  
  // 4th var
  varmap["cp0"]=4;
  cp0=MASA_VAR_DEFAULT;
  vararr.push_back(&cp0);   

  // 5th var
  varmap["at"]=5;
  at=MASA_VAR_DEFAULT;
  vararr.push_back(&at);

  // 6th var
  varmap["rho"]=6;
  rho=MASA_VAR_DEFAULT;
  vararr.push_back(&rho);

  // 7th var
  varmap["by"]=7;
  by=MASA_VAR_DEFAULT;
  vararr.push_back(&by);

  // 8th var
  varmap["bt"]=7;
  bt=MASA_VAR_DEFAULT;
  vararr.push_back(&bt);

  // 7th var
  varmap["cz"]=8;
  cz=MASA_VAR_DEFAULT;
  vararr.push_back(&cz);

  // 8th var
  varmap["ct"]=9;
  ct=MASA_VAR_DEFAULT;
  vararr.push_back(&ct);

}//done with constructor


MASA::heateq_1d_unsteady_var::heateq_1d_unsteady_var()
{
  mmsname = "heateq_1d_unsteady_var";
  dimension=1;

  //first variable (dummy) "axp" -- load map and array
  axp=MASA_VAR_DEFAULT;
  vararr.push_back(&axp);

  // initalize other variables
  varmap["ax"]=1;                   // incriment map location
  ax=MASA_VAR_DEFAULT;              // need to initialize all variables!
  vararr.push_back(&ax);            // add variable to pointer array
  
  // 2nd var
  varmap["k0"]=2;
  k0=MASA_VAR_DEFAULT;
  vararr.push_back(&k0);   
  
  // 3rd var
  varmap["dt"]=3;
  dt=MASA_VAR_DEFAULT;
  vararr.push_back(&dt);   
  
  // 4th var
  varmap["cp0"]=4;
  cp0=MASA_VAR_DEFAULT;
  vararr.push_back(&cp0);   

  // 5th var
  varmap["at"]=5;
  at=MASA_VAR_DEFAULT;
  vararr.push_back(&at);

  // 6th var
  varmap["rho"]=6;
  rho=MASA_VAR_DEFAULT;
  vararr.push_back(&rho);

  // 7th var
  varmap["cp1"]=7;
  cp1=MASA_VAR_DEFAULT;
  vararr.push_back(&cp1);

  // 8th var
  varmap["cp2"]=8;
  cp2=MASA_VAR_DEFAULT;
  vararr.push_back(&cp2);

  // 9th var
  varmap["k1"]=9;
  k1=MASA_VAR_DEFAULT;
  vararr.push_back(&k1);   

  // 10th var
  varmap["k2"]=10;
  k2=MASA_VAR_DEFAULT;
  vararr.push_back(&k2);   

}//done with constructor


MASA::heateq_2d_unsteady_var::heateq_2d_unsteady_var()
{
  mmsname = "heateq_2d_unsteady_var";
  dimension=2;


  //first variable (dummy) "axp" -- load map and array
  axp=MASA_VAR_DEFAULT;
  vararr.push_back(&axp);

  // initalize other variables
  varmap["ax"]=1;                   // incriment map location
  ax=MASA_VAR_DEFAULT;              // need to initialize all variables!
  vararr.push_back(&ax);            // add variable to pointer array
  
  // 2nd var
  varmap["k0"]=2;
  k0=MASA_VAR_DEFAULT;
  vararr.push_back(&k0);   
  
  // 3rd var
  varmap["dt"]=3;
  dt=MASA_VAR_DEFAULT;
  vararr.push_back(&dt);   
  
  // 4th var
  varmap["cp0"]=4;
  cp0=MASA_VAR_DEFAULT;
  vararr.push_back(&cp0);   

  // 5th var
  varmap["at"]=5;
  at=MASA_VAR_DEFAULT;
  vararr.push_back(&at);

  // 6th var
  varmap["rho"]=6;
  rho=MASA_VAR_DEFAULT;
  vararr.push_back(&rho);

  // 7th var
  varmap["cp1"]=7;
  cp1=MASA_VAR_DEFAULT;
  vararr.push_back(&cp1);

  // 8th var
  varmap["cp2"]=8;
  cp2=MASA_VAR_DEFAULT;
  vararr.push_back(&cp2);

  // 9th var
  varmap["k1"]=9;
  k1=MASA_VAR_DEFAULT;
  vararr.push_back(&k1);   

  // 10th var
  varmap["k2"]=10;
  k2=MASA_VAR_DEFAULT;
  vararr.push_back(&k2);   

  // 11th var
  varmap["by"]=11;
  by=MASA_VAR_DEFAULT;
  vararr.push_back(&by);   

  // 12th var
  varmap["bt"]=12;
  bt=MASA_VAR_DEFAULT;
  vararr.push_back(&bt);

}//done with constructor


MASA::heateq_3d_unsteady_var::heateq_3d_unsteady_var()
{
  mmsname = "heateq_3d_unsteady_var";
  dimension=3;

  //first variable (dummy) "axp" -- load map and array
  axp=MASA_VAR_DEFAULT;
  vararr.push_back(&axp);

  // initalize other variables
  varmap["ax"]=1;                   // incriment map location
  ax=MASA_VAR_DEFAULT;              // need to initialize all variables!
  vararr.push_back(&ax);            // add variable to pointer array
  
  // 2nd var
  varmap["k0"]=2;
  k0=MASA_VAR_DEFAULT;
  vararr.push_back(&k0);   
  
  // 3rd var
  varmap["dt"]=3;
  dt=MASA_VAR_DEFAULT;
  vararr.push_back(&dt);   
  
  // 4th var
  varmap["cp0"]=4;
  cp0=MASA_VAR_DEFAULT;
  vararr.push_back(&cp0);   

  // 5th var
  varmap["at"]=5;
  at=MASA_VAR_DEFAULT;
  vararr.push_back(&at);

  // 6th var
  varmap["rho"]=6;
  rho=MASA_VAR_DEFAULT;
  vararr.push_back(&rho);

  // 7th var
  varmap["cp1"]=7;
  cp1=MASA_VAR_DEFAULT;
  vararr.push_back(&cp1);

  // 8th var
  varmap["cp2"]=8;
  cp2=MASA_VAR_DEFAULT;
  vararr.push_back(&cp2);

  // 9th var
  varmap["k1"]=9;
  k1=MASA_VAR_DEFAULT;
  vararr.push_back(&k1);   

  // 10th var
  varmap["k2"]=10;
  k2=MASA_VAR_DEFAULT;
  vararr.push_back(&k2);   

  // 11th var
  varmap["by"]=11;
  by=MASA_VAR_DEFAULT;
  vararr.push_back(&by);   

  // 12th var
  varmap["bt"]=12;
  bt=MASA_VAR_DEFAULT;
  vararr.push_back(&bt);

  // 13th var
  varmap["cz"]=13;
  cz=MASA_VAR_DEFAULT;
  vararr.push_back(&cz);   

  // 14th var
  varmap["ct"]=14;
  ct=MASA_VAR_DEFAULT;
  vararr.push_back(&ct);

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
  dimension=1;
  
  //first variable (dummy) "axp" -- load map and array
  axp=MASA_VAR_DEFAULT;
  vararr.push_back(&axp);
  
  // initalize other variables
  varmap["ax"]=1;
  ax=MASA_VAR_DEFAULT;              // need to initialize all variables!
  vararr.push_back(&ax);

  // 2nd var
  varmap["k0"]=2;
  k0=MASA_VAR_DEFAULT;
  vararr.push_back(&k0);   

  // 3rd var
  varmap["k1"]=3;
  k1=MASA_VAR_DEFAULT;
  vararr.push_back(&k1);

  // 4th var
  varmap["k2"]=4;
  k2=MASA_VAR_DEFAULT;
  vararr.push_back(&k2);

}//done with constructor


MASA::heateq_2d_steady_var::heateq_2d_steady_var()
{
  mmsname = "heateq_2d_steady_var";
  dimension=2;

  //first variable (dummy) "axp" -- load map and array
  axp=MASA_VAR_DEFAULT;
  vararr.push_back(&axp);
  
  // initalize other variables
  varmap["ax"]=1;
  ax=MASA_VAR_DEFAULT;              // need to initialize all variables!
  vararr.push_back(&ax);

  // 2nd var
  varmap["k0"]=2;
  k0=MASA_VAR_DEFAULT;
  vararr.push_back(&k0);   

  // 3rd var
  varmap["k1"]=3;
  k1=MASA_VAR_DEFAULT;
  vararr.push_back(&k1);

  // 4th var
  varmap["k2"]=4;
  k2=MASA_VAR_DEFAULT;
  vararr.push_back(&k2);

  // 5th var
  varmap["by"]=5;
  by=MASA_VAR_DEFAULT;
  vararr.push_back(&by);

}//done with constructor

MASA::heateq_3d_steady_var::heateq_3d_steady_var()
{
  mmsname = "heateq_3d_steady_var";
  dimension=3;

  //first variable (dummy) "axp" -- load map and array
  axp=MASA_VAR_DEFAULT;
  vararr.push_back(&axp);
  
  // initalize other variables
  varmap["ax"]=1;
  ax=MASA_VAR_DEFAULT;              // need to initialize all variables!
  vararr.push_back(&ax);

  // 2nd var
  varmap["k0"]=2;
  k0=MASA_VAR_DEFAULT;
  vararr.push_back(&k0);   

  // 3rd var
  varmap["k1"]=3;
  k1=MASA_VAR_DEFAULT;
  vararr.push_back(&k1);

  // 4th var
  varmap["k2"]=4;
  k2=MASA_VAR_DEFAULT;
  vararr.push_back(&k2);

  // 5th var
  varmap["by"]=5;
  by=MASA_VAR_DEFAULT;
  vararr.push_back(&by);

  // 6th var
  varmap["cz"]=6;
  cz=MASA_VAR_DEFAULT;
  vararr.push_back(&cz);

}//done with constructor
