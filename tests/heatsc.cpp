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

const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

double SourceQ_t_1d(double x, double A_x, double k_0)
{
  double Q_T = A_x * A_x * k_0 * cos(A_x * x);
  return Q_T;
}

double SourceQ_t_2d (
  double x,
  double y,
  double A_x,
  double B_y,
  double k_0)
{
  double Q_T = k_0 * cos(A_x * x) * cos(B_y * y) * (A_x * A_x + B_y * B_y);
  return Q_T;
}

double SourceQ_t_3d (
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

int main()
{
  double tfield,tfield2;
  double param=1.2;
  double x=.5;
  double y=.4;
  double z=.3;

  /// -----------------------------------------------------------------------
  // initalize 1D
  /// -----------------------------------------------------------------------
  masa_init("temp-test-1d","heateq_1d_steady_const");

  masa_init_param();

  double A_x = masa_get_param("A_x");
  double k_0 = masa_get_param("k_0");

  // evaluate source terms (1D)
  masa_sanity_check();
  tfield    = masa_eval_t_source(x);
  tfield2   = SourceQ_t_1d(x,A_x,k_0);

  tfield=fabs(tfield-tfield2);

  if(tfield > threshold)
    {
      cout << "\nMASA REGRESSION TEST FAILED: Heat Equation 1d Steady Constant\n";
      cout << "T Source Term\n";
      cout << "Exceeded Threshold by: " << tfield << endl;
      exit(1);
    }

  // cout << "1D Steady Constant Heat Equation: PASSED\n";

  /// -----------------------------------------------------------------------
  // initalize 2D
  /// -----------------------------------------------------------------------
  masa_init("temp-test-2d","heateq_2d_steady_const");
  masa_init_param();

  A_x = masa_get_param("A_x");
  k_0 = masa_get_param("k_0");
  double B_y = masa_get_param("B_y");

  // evaluate source terms (1D)
  masa_sanity_check();
  tfield    = masa_eval_t_source(x,y);
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

  masa_init("temp-test-3d","heateq_3d_steady_const");
  masa_init_param();

  A_x = masa_get_param("A_x");
  k_0 = masa_get_param("k_0");
  B_y = masa_get_param("B_y");
  double C_z = masa_get_param("C_z");

  // evaluate source terms (1D)
  masa_sanity_check();
  tfield  = masa_eval_t_source(x,y,z);
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
