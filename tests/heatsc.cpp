// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010 The PECOS Development Team
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
// $Author$
// $Id$
//
// heatsc.cpp: tests steady constant heat equation against known source term
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <math.h>
#include <masa.h>
#include <iostream>
#include <stdlib.h>

using namespace std;
using namespace MASA;

typedef double Scalar;

const Scalar threshold = 1.0e-15; // should be small enough to catch any obvious problems

Scalar SourceQ_t_1d(Scalar x, Scalar A_x, Scalar k_0)
{
  Scalar Q_T = A_x * A_x * k_0 * cos(A_x * x);
  return Q_T;
}

Scalar Source_t_1d_an(Scalar A_x,Scalar x)
{
  Scalar T_an;
  T_an = cos(A_x * x);
  return T_an;
}

Scalar SourceQ_t_2d (
  Scalar x,
  Scalar y,
  Scalar A_x,
  Scalar B_y,
  Scalar k_0)
{
  Scalar Q_T = k_0 * cos(A_x * x) * cos(B_y * y) * (A_x * A_x + B_y * B_y);
  return Q_T;
}

Scalar SourceQ_t_3d (
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar A_x,
  Scalar B_y,
  Scalar C_z,
  Scalar k_0)
{
  Scalar Q_T = k_0 * cos(A_x * x) * cos(B_y * y) * cos(C_z * z) * (A_x * A_x + B_y * B_y + C_z * C_z);
  return Q_T;
}

int main()
{
  Scalar tfield,tfield2,tfield3;
  Scalar t_an,t_an2,t_an3;
  Scalar param=1.2;
  Scalar x=.5;
  Scalar y=.4;
  Scalar z=.3;

  /// -----------------------------------------------------------------------
  // initalize 1D
  /// -----------------------------------------------------------------------
  masa_init<Scalar>("temp-test-1d","heateq_1d_steady_const");
  masa_init_param<Scalar>();

  Scalar A_x = masa_get_param<Scalar>("A_x");
  Scalar k_0 = masa_get_param<Scalar>("k_0");

  // evaluate source terms (1D)
  int err = masa_sanity_check<Scalar>();
  if(err != 0)
    {
      cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }
  

  tfield    = masa_eval_t_source<Scalar>(x);
  t_an      = masa_eval_t_an<Scalar>(x);

  tfield2   = SourceQ_t_1d(x,A_x,k_0);
  t_an2     = Source_t_1d_an(A_x,x);

  tfield3 = fabs(tfield-tfield2);
  t_an3   = fabs(t_an-t_an2);

  if(tfield3 > threshold)
    {
      cout << "\nMASA REGRESSION TEST FAILED: Heat Equation 1d Steady Constant\n";
      cout << "T Source Term\n";
      cout << "Exceeded Threshold by: " << tfield3 << endl;
      exit(1);
    }

  if(t_an3 > threshold)
    {
      cout << "\nMASA REGRESSION TEST FAILED: Heat Equation 1d Steady Constant\n";
      cout << "T Analytical Term\n";
      cout << "Exceeded Threshold by: " << t_an3 << endl;
      exit(1);
    }

  // cout << "1D Steady Constant Heat Equation: PASSED\n";

  /// -----------------------------------------------------------------------
  // initalize 2D
  /// ----------------------------------------------------------------------
  masa_init<Scalar>("temp-test-2d","heateq_2d_steady_const");
  masa_init_param<Scalar>();

  A_x = masa_get_param<Scalar>("A_x");
  k_0 = masa_get_param<Scalar>("k_0");
  Scalar B_y = masa_get_param<Scalar>("B_y");

  // evaluate source terms (1D)
  masa_sanity_check<Scalar>();
  tfield    = masa_eval_t_source<Scalar>(x,y);
  tfield2   = SourceQ_t_2d(x,y,A_x,B_y,k_0);

  tfield=fabs(tfield-tfield2);

  if(tfield > threshold)
    {
      cout << "\nMASA REGRESSION TEST FAILED: Heat Equation 2d Steady Constant\n";
      cout << "T Source Term\n";
      cout << "Exceeded Threshold by: " << tfield << endl;
      exit(1);
    }

  // cout << "2D Steady Constant Heat Equation: PASSED\n";

  /// -----------------------------------------------------------------------
  // initalize 3D
  /// -----------------------------------------------------------------------

  masa_init<Scalar>("temp-test-3d","heateq_3d_steady_const");
  string str1;
  masa_init_param<Scalar>();

  A_x = masa_get_param<Scalar>("A_x");
  k_0 = masa_get_param<Scalar>("k_0");
  B_y = masa_get_param<Scalar>("B_y");  
  Scalar C_z = masa_get_param<Scalar>("C_z");

  // evaluate source terms (1D)
  masa_sanity_check<Scalar>();

  tfield  = masa_eval_t_source<Scalar>(x,y,z);
  tfield2 = SourceQ_t_3d(x,y,z,A_x,B_y,C_z,k_0);
  tfield=fabs(tfield-tfield2);

  if(tfield > threshold)
    {
      cout << "\nMASA REGRESSION TEST FAILED: Heat Equation 3d Steady Constant\n";
      cout << "T Source Term\n";
      cout << "Exceeded Threshold by: " << tfield << endl;
      exit(1);
    }

  //cout << "3D Steady Constant Heat Equation: PASSED\n";
  // presumably, all tests passed
  return 0;

}
