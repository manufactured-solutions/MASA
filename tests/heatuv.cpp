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
// heatuv.cpp: test unsteady variable heat equation against known source term
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

double SourceQ_t_1d (
  double x,
  double t,
  double A_x,
  double A_t,
  double D_t,
  double rho,
  double k_0,
  double k_1,
  double k_2,
  double cp_0,
  double cp_1,
  double cp_2)
{
  double Q_T;
  double T_an;
  double *gradT_an;
  Q_T = -pow(cos(D_t * t), 0.2e1) * A_x * A_x * k_1 + 0.3e1 * pow(cos(A_x * x + A_t * t), 0.3e1) * pow(cos(D_t * t), 0.3e1) * A_x * A_x * k_2 - sin(A_x * x + A_t * t) * cos(D_t * t) * A_t * rho * cp_0 - cos(A_x * x + A_t * t) * sin(D_t * t) * D_t * rho * cp_0 + (A_x * A_x * k_0 - 0.2e1 * pow(cos(D_t * t), 0.2e1) * A_x * A_x * k_2 - sin(A_x * x + A_t * t) * cos(D_t * t) * A_t * rho * cp_1 - cos(A_x * x + A_t * t) * sin(D_t * t) * D_t * rho * cp_1) * cos(A_x * x + A_t * t) * cos(D_t * t) + (0.2e1 * A_x * A_x * k_1 - sin(A_x * x + A_t * t) * cos(D_t * t) * A_t * rho * cp_2 - cos(A_x * x + A_t * t) * sin(D_t * t) * D_t * rho * cp_2) * pow(cos(A_x * x + A_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1);
  // gradT_an[0] = -A_x * cos(D_t * t) * sin(A_x * x + A_t * t);
  return Q_T;
}

double SourceQ_t_2d (
  double x,
  double y,
  double t,
  double A_x,
  double A_t,
  double B_y,
  double B_t,
  double D_t,
  double rho,
  double k_0,
  double k_1,
  double k_2,
  double cp_0,
  double cp_1,
  double cp_2)
{
  double Q_T;
  double T_an;
  double *gradT_an;
  Q_T = -pow(cos(B_y * y + B_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1) * A_x * A_x * k_1 - sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(D_t * t) * A_t * rho * cp_0 - pow(cos(A_x * x + A_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1) * B_y * B_y * k_1 - cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(D_t * t) * B_t * rho * cp_0 - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(D_t * t) * D_t * rho * cp_0 + 0.3e1 * pow(cos(A_x * x + A_t * t), 0.3e1) * pow(cos(B_y * y + B_t * t), 0.3e1) * pow(cos(D_t * t), 0.3e1) * (A_x * A_x + B_y * B_y) * k_2 + (0.2e1 * A_x * A_x * k_1 - sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(D_t * t) * A_t * rho * cp_2 + 0.2e1 * B_y * B_y * k_1 - cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(D_t * t) * B_t * rho * cp_2 - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(D_t * t) * D_t * rho * cp_2) * pow(cos(A_x * x + A_t * t), 0.2e1) * pow(cos(B_y * y + B_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1) + (A_x * A_x * k_0 - 0.2e1 * pow(cos(B_y * y + B_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1) * A_x * A_x * k_2 - sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(D_t * t) * A_t * rho * cp_1 + B_y * B_y * k_0 - 0.2e1 * pow(cos(A_x * x + A_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1) * B_y * B_y * k_2 - cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(D_t * t) * B_t * rho * cp_1 - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(D_t * t) * D_t * rho * cp_1) * cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(D_t * t);
  T_an = cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(D_t * t);
  // gradT_an[0] = -A_x * cos(D_t * t) * cos(B_y * y + B_t * t) * sin(A_x * x + A_t * t);
  // gradT_an[1] = -B_y * cos(D_t * t) * cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t);
  return Q_T;
}

double SourceQ_t_3d (
  double x,
  double y,
  double z,
  double t,
  double A_x,
  double A_t,
  double B_y,
  double B_t,
  double C_z,
  double C_t,
  double D_t,
  double rho,
  double k_0,
  double k_1,
  double k_2,
  double cp_0,
  double cp_1,
  double cp_2)
{
  double Q_T;
  double T_an;
  double *gradT_an;
  Q_T = -sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_0 * A_t - cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_0 * B_t - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_0 * C_t - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * sin(D_t * t) * rho * cp_0 * D_t - pow(cos(B_y * y + B_t * t), 0.2e1) * pow(cos(C_z * z + C_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1) * k_1 * A_x * A_x - pow(cos(A_x * x + A_t * t), 0.2e1) * pow(cos(C_z * z + C_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1) * k_1 * B_y * B_y - pow(cos(A_x * x + A_t * t), 0.2e1) * pow(cos(B_y * y + B_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1) * k_1 * C_z * C_z + 0.3e1 * pow(cos(A_x * x + A_t * t), 0.3e1) * pow(cos(B_y * y + B_t * t), 0.3e1) * pow(cos(C_z * z + C_t * t), 0.3e1) * pow(cos(D_t * t), 0.3e1) * (A_x * A_x + B_y * B_y + C_z * C_z) * k_2 + (-sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_1 * A_t - cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_1 * B_t - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_1 * C_t - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * sin(D_t * t) * rho * cp_1 * D_t + k_0 * A_x * A_x + k_0 * B_y * B_y + k_0 * C_z * C_z - 0.2e1 * pow(cos(B_y * y + B_t * t), 0.2e1) * pow(cos(C_z * z + C_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1) * k_2 * A_x * A_x - 0.2e1 * pow(cos(A_x * x + A_t * t), 0.2e1) * pow(cos(C_z * z + C_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1) * k_2 * B_y * B_y - 0.2e1 * pow(cos(A_x * x + A_t * t), 0.2e1) * pow(cos(B_y * y + B_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1) * k_2 * C_z * C_z) * cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) + (-sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_2 * A_t - cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_2 * B_t - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_2 * C_t - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * sin(D_t * t) * rho * cp_2 * D_t + 0.2e1 * k_1 * A_x * A_x + 0.2e1 * k_1 * B_y * B_y + 0.2e1 * k_1 * C_z * C_z) * pow(cos(A_x * x + A_t * t), 0.2e1) * pow(cos(B_y * y + B_t * t), 0.2e1) * pow(cos(C_z * z + C_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1);
  T_an = cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t);
  // gradT_an[0] = -A_x * cos(D_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * sin(A_x * x + A_t * t);
  // gradT_an[1] = -B_y * cos(D_t * t) * cos(A_x * x + A_t * t) * cos(C_z * z + C_t * t) * sin(B_y * y + B_t * t);
  // gradT_an[2] = -C_z * cos(D_t * t) * cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(C_z * z + C_t * t);
  return Q_T;
}

int main()
{
  double tfield,tfield2;
  double param=1.2;
  double t= 1;
  double x=.5;
  double y=.4;
  double z=.3;

  /// -----------------------------------------------------------------------
  // initalize 1D
  /// -----------------------------------------------------------------------
  masa_init("temp-test-1d","heateq_1d_unsteady_var");
  
  masa_init_param();

  masa_set_param("A_x",param);
  double A_x=param;

  masa_set_param("A_t",param);
  double A_t=param;

  masa_set_param("D_t",param);
  double D_t=param;

  masa_set_param("cp_0",param);
  double cp_0=param;

  masa_set_param("cp_1",param);
  double cp_1=param;

  masa_set_param("cp_2",param);
  double cp_2=param;

  masa_set_param("rho",param);
  double rho=param;

  masa_set_param("k_0",param);
  double k_0=param;

  masa_set_param("k_1",param);
  double k_1=param;

  masa_set_param("k_2",param);
  double k_2=param;

  // evaluate source terms (1D)
  int err = masa_sanity_check();
  if(err != 0)
    {
      cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }
  
  tfield   = masa_eval_t_source(x,t);
  tfield2  = SourceQ_t_1d(x,t,A_x,A_t,D_t,rho,k_0,k_1,k_2,cp_0,cp_1,cp_2);

  tfield=fabs(tfield-tfield2);

  if(tfield > threshold)
    {
      cout << "\nMASA REGRESSION TEST FAILED: Heat Equation 1d Unsteady Variable Coefficients\n";
      cout << "T Source Term\n";
      cout << "Exceeded Threshold by: " << tfield << endl;
      exit(1);
    }

  //cout << "1D Unsteady Variable Coefficients Heat Equation: PASSED\n";

  /// -----------------------------------------------------------------------
  // initalize 2D
  /// -----------------------------------------------------------------------
  masa_init("temp-test-2d","heateq_2d_unsteady_var");
  masa_init_param();

  masa_set_param("A_x",param);
  A_x=param;

  masa_set_param("A_t",param);
  A_t=param;

  masa_set_param("D_t",param);
  D_t=param;

  masa_set_param("cp_0",param);
  cp_0=param;

  masa_set_param("cp_1",param);
  cp_1=param;

  masa_set_param("cp_2",param);
  cp_2=param;

  masa_set_param("rho",param);
  rho=param;

  masa_set_param("k_0",param);
  k_0=param;

  masa_set_param("k_1",param);
  k_1=param;

  masa_set_param("k_2",param);
  k_2=param;

  masa_set_param("B_y",param);
  double B_y=param;

  masa_set_param("B_t",param);
  double B_t=param;

  // evaluate source terms (2D)
  err = masa_sanity_check();
  if(err != 0)
    {
      cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }
  
  tfield   = masa_eval_t_source(x,y,t);
  tfield2   = SourceQ_t_2d(x,y,t,A_x,A_t,B_y,B_t,D_t,rho,k_0,k_1,k_2,cp_0,cp_1,cp_2);

  tfield=fabs(tfield-tfield2);

  if(tfield > threshold)
    {
      cout << "\nMASA REGRESSION TEST FAILED: Heat Equation 2d Unsteady Variable Coefficients\n";
      cout << "T Source Term\n";
      cout << "Exceeded Threshold by: " << tfield << endl;
      exit(1);
    }

  //cout << "2D Unsteady Variable Coefficients Heat Equation: PASSED\n";


  /// -----------------------------------------------------------------------
  // initalize 3D
  /// -----------------------------------------------------------------------
  masa_init("temp-test-3d","heateq_3d_unsteady_var");
  masa_init_param();

  masa_set_param("A_x",param);
  A_x=param;

  masa_set_param("A_t",param);
  A_t=param;

  masa_set_param("D_t",param);
  D_t=param;

  masa_set_param("cp_0",param);
  cp_0=param;

  masa_set_param("rho",param);
  rho=param;

  masa_set_param("k_0",param);
  k_0=param;

  masa_set_param("B_y",param);
  B_y=param;

  masa_set_param("B_t",param);
  B_t=param;

  masa_set_param("k_1",param);
  k_1=param;

  masa_set_param("k_2",param);
  k_2=param;

  masa_set_param("cp_0",param);
  cp_0=param;

  masa_set_param("cp_1",param);
  cp_1=param;

  masa_set_param("cp_2",param);
  cp_2=param;

  masa_set_param("C_z",param);
  double C_z=param;

  masa_set_param("C_t",param);
  double C_t=param;

  // evaluate source terms (3D)
  err = masa_sanity_check();
  if(err != 0)
    {
      cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }
  
  tfield   = masa_eval_t_source(x,y,z,t);
  tfield2   = SourceQ_t_3d(x,y,z,t,A_x,A_t,B_y,B_t,C_z,C_t,D_t,rho,k_0,k_1,k_2,cp_0,cp_1,cp_2);

  tfield=fabs(tfield-tfield2);

  if(tfield > threshold)
    {
      cout << "\nMASA REGRESSION TEST FAILED: Heat Equation 3d Unsteady Variable Coefficients\n";
      cout << "T Source Term\n";
      cout << "Exceeded Threshold by: " << tfield << endl;
      exit(1);
    }

  //cout << "3D Unsteady Variable Coefficients Heat Equation: PASSED\n";

  // presumably, all tests passed
  return 0;

}
