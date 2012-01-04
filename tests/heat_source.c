// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012 The PECOS Development Team
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
// $Id: euler-source.c 12693 2010-08-26 03:35:34Z nick $
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <math.h>
#include <stdio.h>

double eval_1d_exact_t_(double x,double A_x)
{  
  double exact_t;
  exact_t = cos(A_x * x);
  return exact_t;
}
 
double eval_1d_source_t_(double x, double A_x, double k_0)
{
  double Q_T = A_x * A_x * k_0 * cos(A_x * x);
  return Q_T;
}

//--------------------------------------------------------------------------
// 
// 
//    2D problem
// 
//
//--------------------------------------------------------------------------

  
double eval_2d_exact_t_(double x,double y,double A_x,double B_y)
{
  double exact_t;
  exact_t = cos(A_x * x) * cos(B_y * y);
  return exact_t;
}

double eval_2d_source_t_(
  double x,
  double y,
  double A_x,
  double B_y,
  double k_0)
{
  double Q_T = k_0 * cos(A_x * x) * cos(B_y * y) * (A_x * A_x + B_y * B_y);
  return Q_T;
}

//--------------------------------------------------------------------------
// 
// 
//    3D problem
// 
//
//--------------------------------------------------------------------------
  
double eval_3d_exact_t_(double x,double y,double z,double A_x,double B_y,double C_z)
{
  double exact_t;
  exact_t = cos(A_x * x) * cos(B_y * y) * cos(C_z * z);
  return exact_t;
}

double eval_3d_source_t_(
  double x,
  double y,
  double z,
  double A_x,
  double B_y,
  double C_z,
  double k_0)
{
  double Q_T = k_0 * cos(A_x * x) * cos(B_y * y) * cos(C_z * z) * (A_x * A_x + B_y * B_y + C_z * C_z);
  return Q_T;
}
