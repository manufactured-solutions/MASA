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

    register_var("A_x",&A_x);   
    register_var("k_0",&k_0);

}//done with constructor

double MASA::heateq_1d_steady_const::eval_q_t(double x)
{
  double Q_T;
  Q_T = A_x * A_x * k_0 * cos(A_x * x);
  return Q_T;
}

double MASA::heateq_1d_steady_const::eval_an(double x)
{
  double T_an;
  T_an = cos(A_x * x);
  return T_an;
}

MASA::heateq_2d_steady_const::heateq_2d_steady_const()
{
    mmsname = "heateq_2d_steady_const";
    dimension=2;

    register_var("A_x",&A_x);   
    register_var("k_0",&k_0);
    register_var("B_y",&B_y);
    
}//done with constructor

double MASA::heateq_2d_steady_const::eval_q_t(double x,double y)
{
  double Q_T = k_0 * cos(A_x * x) * cos(B_y * y) * (A_x * A_x + B_y * B_y);
  return Q_T;
}


MASA::heateq_3d_steady_const::heateq_3d_steady_const()
{
    mmsname = "heateq_3d_steady_const";
    dimension=3;

    register_var("A_x",&A_x);   
    register_var("k_0",&k_0);
    register_var("B_y",&B_y);
    register_var("C_z",&C_z);

}//done with constructor

double MASA::heateq_3d_steady_const::eval_q_t(double x,double y,double z)
{
  double Q_T = k_0 * cos(A_x * x) * cos(B_y * y) * cos(C_z * z) * (A_x * A_x + B_y * B_y + C_z * C_z);
  return Q_T;
}



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

  register_var("ax",&ax);   
  register_var("k0",&k0);
  register_var("dt",&dt);
  register_var("cp0",&cp0);
  register_var("at",&at);
  register_var("rho",&rho);

}//done with constructor


MASA::heateq_2d_unsteady_const::heateq_2d_unsteady_const()
{
  mmsname = "heateq_2d_unsteady_const";
  dimension=2;
    
  register_var("ax",&ax);   
  register_var("k0",&k0);
  register_var("dt",&dt);
  register_var("cp0",&cp0);
  register_var("at",&at);
  register_var("rho",&rho);
  register_var("by",&by);
  register_var("bt",&bt);

}//done with constructor


MASA::heateq_3d_unsteady_const::heateq_3d_unsteady_const()
{
  mmsname = "heateq_3d_unsteady_const";
  dimension=3;
  
  register_var("ax",&ax);   
  register_var("k0",&k0);
  register_var("dt",&dt);
  register_var("cp0",&cp0);
  register_var("at",&at);
  register_var("rho",&rho);
  register_var("by",&by);
  register_var("bt",&bt);
  register_var("cz",&cz);
  register_var("ct",&ct);

}//done with constructor


MASA::heateq_1d_unsteady_var::heateq_1d_unsteady_var()
{
  mmsname = "heateq_1d_unsteady_var";
  dimension=1;

  register_var("ax",&ax);   
  register_var("k0",&k0);
  register_var("dt",&dt);
  register_var("cp0",&cp0);
  register_var("at",&at);
  register_var("rho",&rho);
  register_var("cp1",&cp1);
  register_var("cp2",&cp2);
  register_var("k1",&k1);
  register_var("k2",&k2);

}//done with constructor


MASA::heateq_2d_unsteady_var::heateq_2d_unsteady_var()
{
  mmsname = "heateq_2d_unsteady_var";
  dimension=2;

  register_var("ax",&ax);   
  register_var("k0",&k0);
  register_var("dt",&dt);
  register_var("cp0",&cp0);
  register_var("at",&at);
  register_var("rho",&rho);
  register_var("cp1",&cp1);
  register_var("cp2",&cp2);
  register_var("k1",&k1);
  register_var("k2",&k2);
  register_var("by",&by);
  register_var("bt",&bt);

}//done with constructor


MASA::heateq_3d_unsteady_var::heateq_3d_unsteady_var()
{
  mmsname = "heateq_3d_unsteady_var";
  dimension=3;

  register_var("ax",&ax);   
  register_var("k0",&k0);
  register_var("dt",&dt);
  register_var("cp0",&cp0);
  register_var("at",&at);
  register_var("rho",&rho);
  register_var("cp1",&cp1);
  register_var("cp2",&cp2);
  register_var("k1",&k1);
  register_var("k2",&k2);
  register_var("by",&by);
  register_var("bt",&bt);
  register_var("cz",&cz);
  register_var("ct",&ct);

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

  register_var("ax",&ax);   
  register_var("k0",&k0);
  register_var("k1",&k1);
  register_var("k2",&k2);

}//done with constructor


MASA::heateq_2d_steady_var::heateq_2d_steady_var()
{
  mmsname = "heateq_2d_steady_var";
  dimension=2;

  register_var("ax",&ax);   
  register_var("k0",&k0);
  register_var("k1",&k1);
  register_var("k2",&k2);
  register_var("by",&by);

}//done with constructor

MASA::heateq_3d_steady_var::heateq_3d_steady_var()
{
  mmsname = "heateq_3d_steady_var";
  dimension=3;

  register_var("ax",&ax);   
  register_var("k0",&k0);
  register_var("k1",&k1);
  register_var("k2",&k2);
  register_var("by",&by);
  register_var("cz",&cz);

}//done with constructor
