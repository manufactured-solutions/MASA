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

  register_var("A_x",&A_x);   
  register_var("k_0",&k_0);
  register_var("D_t",&D_t);
  register_var("cp_0",&cp_0);
  register_var("A_t",&A_t);
  register_var("rho",&rho);

}//done with constructor

double MASA::heateq_1d_unsteady_const::eval_q_t(double x,double t)
{
  double Q_T = cos(A_x * x + A_t * t) * cos(D_t * t) * k_0 * A_x * A_x - (sin(A_x * x + A_t * t) * cos(D_t * t) * A_t + cos(A_x * x + A_t * t) * sin(D_t * t) * D_t) * rho * cp_0;
  return Q_T;
}



MASA::heateq_2d_unsteady_const::heateq_2d_unsteady_const()
{
  mmsname = "heateq_2d_unsteady_const";
  dimension=2;
    
  register_var("A_x",&A_x);   
  register_var("k_0",&k_0);
  register_var("D_t",&D_t);
  register_var("cp_0",&cp_0);
  register_var("A_t",&A_t);
  register_var("rho",&rho);
  register_var("B_y",&B_y);
  register_var("B_t",&B_t);

}//done with constructor

double MASA::heateq_2d_unsteady_const::eval_q_t(double x,double y, double t)
{
  double Q_T = -(sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(D_t * t) * A_t + cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(D_t * t) * B_t + cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(D_t * t) * D_t) * rho * cp_0 + (A_x * A_x + B_y * B_y) * cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(D_t * t) * k_0;
  return Q_T;
}


MASA::heateq_3d_unsteady_const::heateq_3d_unsteady_const()
{
  mmsname = "heateq_3d_unsteady_const";
  dimension=3;
  
  register_var("A_x",&A_x);   
  register_var("k_0",&k_0);
  register_var("D_t",&D_t);
  register_var("cp_0",&cp_0);
  register_var("A_t",&A_t);
  register_var("rho",&rho);
  register_var("B_y",&B_y);
  register_var("B_t",&B_t);
  register_var("C_z",&C_z);
  register_var("C_t",&C_t);

}//done with constructor

double MASA::heateq_3d_unsteady_const::eval_q_t(double x,double y, double z,double t)
{
  double Q_T = -(sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * A_t + cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * B_t + cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(C_z * z + C_t * t) * cos(D_t * t) * C_t + cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * sin(D_t * t) * D_t) * rho * cp_0 + (A_x * A_x + B_y * B_y + C_z * C_z) * cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * k_0;
  return Q_T;
}


MASA::heateq_1d_unsteady_var::heateq_1d_unsteady_var()
{
  mmsname = "heateq_1d_unsteady_var";
  dimension=1;

  register_var("A_x",&A_x);   
  register_var("k_0",&k_0);
  register_var("D_t",&D_t);
  register_var("cp_0",&cp_0);
  register_var("A_t",&A_t);
  register_var("rho",&rho);
  register_var("cp_1",&cp_1);
  register_var("cp_2",&cp_2);
  register_var("k_1",&k_1);
  register_var("k_2",&k_2);

}//done with constructor


MASA::heateq_2d_unsteady_var::heateq_2d_unsteady_var()
{
  mmsname = "heateq_2d_unsteady_var";
  dimension=2;

  register_var("A_x",&A_x);   
  register_var("k_0",&k_0);
  register_var("D_t",&D_t);
  register_var("cp_0",&cp_0);
  register_var("A_t",&A_t);
  register_var("rho",&rho);
  register_var("cp_1",&cp_1);
  register_var("cp_2",&cp_2);
  register_var("k_1",&k_1);
  register_var("k_2",&k_2);
  register_var("B_y",&B_y);
  register_var("B_t",&B_t);

}//done with constructor


MASA::heateq_3d_unsteady_var::heateq_3d_unsteady_var()
{
  mmsname = "heateq_3d_unsteady_var";
  dimension=3;

  register_var("A_x",&A_x);   
  register_var("k_0",&k_0);
  register_var("D_t",&D_t);
  register_var("cp_0",&cp_0);
  register_var("A_t",&A_t);
  register_var("rho",&rho);
  register_var("cp_1",&cp_1);
  register_var("cp_2",&cp_2);
  register_var("k_1",&k_1);
  register_var("k_2",&k_2);
  register_var("B_y",&B_y);
  register_var("B_t",&B_t);
  register_var("C_z",&C_z);
  register_var("C_t",&C_t);

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

  register_var("A_x",&A_x);   
  register_var("k_0",&k_0);
  register_var("k_1",&k_1);
  register_var("k_2",&k_2);

}//done with constructor


MASA::heateq_2d_steady_var::heateq_2d_steady_var()
{
  mmsname = "heat_eq_2d_steady_var";
  dimension=2;

  register_var("A_x",&A_x);   
  register_var("k_0",&k_0);
  register_var("k_1",&k_1);
  register_var("k_2",&k_2);
  register_var("B_y",&B_y);

}//done with constructor

MASA::heateq_3d_steady_var::heateq_3d_steady_var()
{
  mmsname = "heat_eq_3d_steady_var";
  dimension=3;

  register_var("A_x",&A_x);   
  register_var("k_0",&k_0);
  register_var("k1",&k_1);
  register_var("k2",&k_2);
  register_var("B_y",&B_y);
  register_var("C_z",&C_z);

}//done with constructor
