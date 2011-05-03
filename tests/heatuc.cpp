// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011 The PECOS Development Team
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
// heatuc.cpp: test unsteady constant heat equation against known source term
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <tests.h>
#include <cmath>

using namespace std;
using namespace MASA;

template<typename Scalar>
Scalar SourceQ_t_1d(
  Scalar x,
  Scalar t,
  Scalar A_x,
  Scalar A_t,
  Scalar D_t,
  Scalar k_0,
  Scalar cp_0,
  Scalar rho)
{
  Scalar Q_T = cos(A_x * x + A_t * t) * cos(D_t * t) * k_0 * A_x * A_x - (sin(A_x * x + A_t * t) * cos(D_t * t) * A_t + cos(A_x * x + A_t * t) * sin(D_t * t) * D_t) * rho * cp_0;
  return Q_T;
}

template<typename Scalar>
Scalar SourceQ_t_2d (
  Scalar x,
  Scalar y,
  Scalar t,
  Scalar A_x,
  Scalar A_t,
  Scalar B_y,
  Scalar B_t,
  Scalar D_t,
  Scalar rho,
  Scalar k_0,
  Scalar cp_0)
{
  Scalar Q_T = -(sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(D_t * t) * A_t + cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(D_t * t) * B_t + cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(D_t * t) * D_t) * rho * cp_0 + (A_x * A_x + B_y * B_y) * cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(D_t * t) * k_0;
  return Q_T;
}

template<typename Scalar>
Scalar SourceQ_t_3d (
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar t,
  Scalar A_x,
  Scalar A_t,
  Scalar B_y,
  Scalar B_t,
  Scalar C_z,
  Scalar C_t,
  Scalar D_t,
  Scalar k_0,
  Scalar cp_0,
  Scalar rho)
{
  Scalar Q_T = -(sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * A_t + cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * B_t + cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(C_z * z + C_t * t) * cos(D_t * t) * C_t + cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * sin(D_t * t) * D_t) * rho * cp_0 + (A_x * A_x + B_y * B_y + C_z * C_z) * cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * k_0;
  return Q_T;
}

template<typename Scalar>
int run_regression()
{
  const Scalar threshold = 5 * numeric_limits<Scalar>::epsilon();

  int err=1;
  Scalar tfield,tfield2;
  Scalar param=1.2;
  Scalar t= 1;
  Scalar x=.5;
  Scalar y=.4;
  Scalar z=.3;

  /// -----------------------------------------------------------------------
  // initalize 1D
  /// -----------------------------------------------------------------------
  masa_init<Scalar>("temp-test-1d","heateq_1d_unsteady_const");
  
  masa_init_param<Scalar>();

  masa_set_param<Scalar>("A_x",param);
  Scalar A_x=param;

  masa_set_param<Scalar>("A_t",param);
  Scalar A_t=param;

  masa_set_param<Scalar>("D_t",param);
  Scalar D_t=param;

  masa_set_param<Scalar>("cp_0",param);
  Scalar cp_0=param;

  masa_set_param<Scalar>("rho",param);
  Scalar rho=param;

  masa_set_param<Scalar>("k_0",param);
  Scalar k_0=param;

  // evaluate source terms (1D)
  err = masa_sanity_check<Scalar>();
  if(err != 0)
    {
      cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }
  
  tfield   = masa_eval_source_t<Scalar>(x,t);
  tfield2  = SourceQ_t_1d(x,t,A_x,A_t,D_t,k_0,cp_0,rho);

  tfield=std::abs(tfield-tfield2);

  threshcheck(tfield,threshold);
  //cout << "1D Unsteady Constant Heat Equation: PASSED\n";

  /// -----------------------------------------------------------------------
  // initalize 2D
  /// -----------------------------------------------------------------------
  masa_init<Scalar>("temp-test-2d","heateq_2d_unsteady_const");

  masa_init_param<Scalar>();

  masa_set_param<Scalar>("A_x",param);
  A_x=param;

  masa_set_param<Scalar>("A_t",param);
  A_t=param;

  masa_set_param<Scalar>("D_t",param);
  D_t=param;

  masa_set_param<Scalar>("cp_0",param);
  cp_0=param;

  masa_set_param<Scalar>("rho",param);
  rho=param;

  masa_set_param<Scalar>("k_0",param);
  k_0=param;

  masa_set_param<Scalar>("B_y",param);
  Scalar B_y=param;

  masa_set_param<Scalar>("B_t",param);
  Scalar B_t=param;

  // evaluate source terms (2D)
  err = masa_sanity_check<Scalar>();
  if(err != 0)
    {
      cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }
  
  tfield   = masa_eval_source_t<Scalar>(x,y,t);
  tfield2   = SourceQ_t_2d(x,y,t,A_x,A_t,B_y,B_t,D_t,rho,k_0,cp_0);

  tfield=std::abs(tfield-tfield2);

  threshcheck(tfield,threshold);

  /// -----------------------------------------------------------------------
  // initalize 3D
  /// -----------------------------------------------------------------------
  masa_init<Scalar>("temp-test-3d","heateq_3d_unsteady_const");

  masa_init_param<Scalar>();

  masa_set_param<Scalar>("A_x",param);
  A_x=param;

  masa_set_param<Scalar>("A_t",param);
  A_t=param;

  masa_set_param<Scalar>("D_t",param);
  D_t=param;

  masa_set_param<Scalar>("cp_0",param);
  cp_0=param;

  masa_set_param<Scalar>("rho",param);
  rho=param;

  masa_set_param<Scalar>("k_0",param);
  k_0=param;

  masa_set_param<Scalar>("B_y",param);
  B_y=param;

  masa_set_param<Scalar>("B_t",param);
  B_t=param;

  masa_set_param<Scalar>("C_z",param);
  Scalar C_z=param;

  masa_set_param<Scalar>("C_t",param);
  Scalar C_t=param;

  // evaluate source terms (3D)
  err = masa_sanity_check<Scalar>();
  if(err != 0)
    {
      cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }
  
  tfield    = masa_eval_source_t<Scalar>(x,y,z,t);
  tfield2   = SourceQ_t_3d(x,y,z,t,A_x,A_t,B_y,B_t,C_z,C_t,D_t,k_0,cp_0,rho);

  tfield=std::abs(tfield-tfield2);

  threshcheck(tfield,threshold);
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
