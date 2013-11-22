// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012,2013 The PECOS Development Team
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
// heatsv.cpp: tests steady variable heat equation against known source term
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <tests.h>
#include <cmath>

using namespace MASA;
using namespace std;

template<typename Scalar>
Scalar SourceQ_t_1d(
  Scalar x,
  Scalar A_x,
  Scalar k_0,
  Scalar k_1,
  Scalar k_2)
{
  Scalar Q_T = Scalar(0.3e1) * A_x * A_x * k_2 * std::pow(cos(A_x * x), Scalar(0.3e1)) + Scalar(0.2e1) * A_x * A_x * k_1 * std::pow(cos(A_x * x), Scalar(0.2e1)) - A_x * A_x * k_1 + (k_0 - Scalar(0.2e1) * k_2) * A_x * A_x * std::cos(A_x * x);
  return Q_T;
}

template<typename Scalar>
Scalar SourceQ_t_2d (
  Scalar x,
  Scalar y,
  Scalar A_x,
  Scalar B_y,
  Scalar k_0,
  Scalar k_1,
  Scalar k_2)
{
  Scalar Q_T = (Scalar(0.3e1) * A_x * A_x + Scalar(0.3e1) * B_y * B_y) * k_2 * std::pow(cos(A_x * x), Scalar(0.3e1)) * std::pow(cos(B_y * y), Scalar(0.3e1)) + (Scalar(0.2e1) * A_x * A_x + Scalar(0.2e1) * B_y * B_y) * k_1 * std::pow(cos(A_x * x), Scalar(0.2e1)) * std::pow(cos(B_y * y), Scalar(0.2e1)) - (std::pow(cos(B_y * y), Scalar(0.2e1)) * A_x * A_x + std::pow(cos(A_x * x), Scalar(0.2e1)) * B_y * B_y) * k_1 + (k_0 * A_x * A_x + k_0 * B_y * B_y - Scalar(0.2e1) * std::pow(cos(B_y * y), Scalar(0.2e1)) * k_2 * A_x * A_x - Scalar(0.2e1) * std::pow(cos(A_x * x), Scalar(0.2e1)) * k_2 * B_y * B_y) * std::cos(A_x * x) * std::cos(B_y * y);
  return Q_T;
}

template<typename Scalar>
Scalar SourceQ_t_3d (
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar A_x,
  Scalar B_y,
  Scalar C_z,
  Scalar k_0,
  Scalar k_1,
  Scalar k_2)
{
  Scalar Q_T = (Scalar(0.3e1) * A_x * A_x + Scalar(0.3e1) * B_y * B_y + Scalar(0.3e1) * C_z * C_z) * k_2 * std::pow(cos(A_x * x), Scalar(0.3e1)) * std::pow(cos(B_y * y), Scalar(0.3e1)) * std::pow(cos(C_z * z), Scalar(0.3e1)) + (Scalar(0.2e1) * A_x * A_x + Scalar(0.2e1) * B_y * B_y + Scalar(0.2e1) * C_z * C_z) * k_1 * std::pow(cos(A_x * x), Scalar(0.2e1)) * std::pow(cos(B_y * y), Scalar(0.2e1)) * std::pow(cos(C_z * z), Scalar(0.2e1)) - (std::pow(cos(B_y * y), Scalar(0.2e1)) * std::pow(cos(C_z * z), Scalar(0.2e1)) * A_x * A_x + std::pow(cos(A_x * x), Scalar(0.2e1)) * std::pow(cos(C_z * z), Scalar(0.2e1)) * B_y * B_y + std::pow(cos(A_x * x), Scalar(0.2e1)) * std::pow(cos(B_y * y), Scalar(0.2e1)) * C_z * C_z) * k_1 + (k_0 * A_x * A_x + k_0 * B_y * B_y + k_0 * C_z * C_z - Scalar(0.2e1) * std::pow(cos(B_y * y), Scalar(0.2e1)) * std::pow(cos(C_z * z), Scalar(0.2e1)) * k_2 * A_x * A_x - Scalar(0.2e1) * std::pow(cos(A_x * x), Scalar(0.2e1)) * std::pow(cos(C_z * z), Scalar(0.2e1)) * k_2 * B_y * B_y - Scalar(0.2e1) * std::pow(cos(A_x * x), Scalar(0.2e1)) * std::pow(cos(B_y * y), Scalar(0.2e1)) * k_2 * C_z * C_z) * std::cos(A_x * x) * std::cos(B_y * y) * std::cos(C_z * z);
  return Q_T;
}

template<typename Scalar>
int run_regression()
{
  Scalar tfield;
  Scalar tfield2;
  Scalar param=1.2;
  Scalar x=.5;
  Scalar y=.4;
  Scalar z=.3;

  /// -----------------------------------------------------------------------
  // initalize 1D
  /// -----------------------------------------------------------------------
  masa_init<Scalar>("temp-test-1d","heateq_1d_steady_var");
  masa_init_param<Scalar>();

  masa_set_param<Scalar>("A_x",param);
  Scalar A_x=param;

  masa_set_param<Scalar>("k_0",param);
  Scalar k_0=param;

  masa_set_param<Scalar>("k_1",param);
  Scalar k_1=param;

  masa_set_param<Scalar>("k_2",param);
  Scalar k_2=param;

  // evaluate source terms (1D)
  int err = masa_sanity_check<Scalar>();
  if(err != 0)
    {
      cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }
  

  tfield    = masa_eval_source_t<Scalar>(x);  
  tfield2   = SourceQ_t_1d(x,A_x,k_0,k_1,k_2);

  tfield=std::abs(tfield-tfield2);

  threshcheck(tfield);

  //cout << "1D Steady Variable Coefficient Heat Equation: PASSED\n";
  //cout << "Residual: "<< tfield << endl;

  /// -----------------------------------------------------------------------
  // initalize 2D
  /// -----------------------------------------------------------------------
  masa_init<Scalar>("temp-test-2d","heateq_2d_steady_var");
  masa_init_param<Scalar>();

  masa_set_param<Scalar>("A_x",param);
  A_x=param;

  masa_set_param<Scalar>("B_y",param);
  Scalar B_y=param;

  masa_set_param<Scalar>("k_0",param);
  k_0=param;

  masa_set_param<Scalar>("k_1",param);
  k_1=param;

  masa_set_param<Scalar>("k_2",param);
  k_2=param;

  // evaluate source terms (2D)
  err = masa_sanity_check<Scalar>();
  if(err != 0)
    {
      cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }
  
  tfield    = masa_eval_source_t<Scalar>(x,y);  
  tfield2   = SourceQ_t_2d(x,y,A_x,B_y,k_0,k_1,k_2);

  tfield=std::abs(tfield-tfield2);


  threshcheck(tfield);

  //cout << "2D Steady Variable Coefficient Heat Equation: PASSED\n";
  //cout << "Residual: "<< tfield << endl;

  /// -----------------------------------------------------------------------
  // initalize 3D
  /// -----------------------------------------------------------------------
  masa_init<Scalar>("temp-test-3d","heateq_3d_steady_var");
  masa_init_param<Scalar>();

  Scalar C_z;

  A_x = masa_get_param<Scalar>("A_x");
  B_y = masa_get_param<Scalar>("B_y");
  C_z = masa_get_param<Scalar>("C_z");
  k_0 = masa_get_param<Scalar>("k_0");
  k_1 = masa_get_param<Scalar>("k_1");
  k_2 = masa_get_param<Scalar>("k_2");

  // evaluate source terms (2D)
  err = masa_sanity_check<Scalar>();
  if(err != 0)
    {
      cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }
  
  tfield    = masa_eval_source_t<Scalar>(x,y,z);  
  tfield2   = SourceQ_t_3d(x,y,z,A_x,B_y,C_z,k_0,k_1,k_2);

  tfield=std::abs(tfield-tfield2);

  threshcheck(tfield);

  // presumably, all tests passed
  return 0;
}

int main()
{
  int err=0;

  err += run_regression<double>();
  err += run_regression<long double>();

  return err;
}
