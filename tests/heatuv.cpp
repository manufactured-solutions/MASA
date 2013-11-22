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
// heatuv.cpp: test unsteady variable heat equation against known source term
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <tests.h>
#include <cmath>

using namespace std;
using namespace MASA;

typedef double Scalar;

const Scalar threshold = 5 * numeric_limits<Scalar>::epsilon();

template<typename Scalar>
Scalar SourceQ_t_1d (
  Scalar x,
  Scalar t,
  Scalar A_x,
  Scalar A_t,
  Scalar D_t,
  Scalar rho,
  Scalar k_0,
  Scalar k_1,
  Scalar k_2,
  Scalar cp_0,
  Scalar cp_1,
  Scalar cp_2)
{
  Scalar Q_T;
  //Scalar *gradexact_t;
  Q_T = -std::pow(cos(D_t * t), 0.2e1) * A_x * A_x * k_1 + 0.3e1 * std::pow(cos(A_x * x + A_t * t), 0.3e1) * std::pow(cos(D_t * t), 0.3e1) * A_x * A_x * k_2 - std::sin(A_x * x + A_t * t) * std::cos(D_t * t) * A_t * rho * cp_0 - std::cos(A_x * x + A_t * t) * std::sin(D_t * t) * D_t * rho * cp_0 + (A_x * A_x * k_0 - 0.2e1 * std::pow(cos(D_t * t), 0.2e1) * A_x * A_x * k_2 - std::sin(A_x * x + A_t * t) * std::cos(D_t * t) * A_t * rho * cp_1 - std::cos(A_x * x + A_t * t) * std::sin(D_t * t) * D_t * rho * cp_1) * std::cos(A_x * x + A_t * t) * std::cos(D_t * t) + (0.2e1 * A_x * A_x * k_1 - std::sin(A_x * x + A_t * t) * std::cos(D_t * t) * A_t * rho * cp_2 - std::cos(A_x * x + A_t * t) * std::sin(D_t * t) * D_t * rho * cp_2) * std::pow(cos(A_x * x + A_t * t), 0.2e1) * std::pow(cos(D_t * t), 0.2e1);
  // gradexact_t[0] = -A_x * std::cos(D_t * t) * std::sin(A_x * x + A_t * t);
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
  Scalar k_1,
  Scalar k_2,
  Scalar cp_0,
  Scalar cp_1,
  Scalar cp_2)
{
  Scalar Q_T;
  //Scalar *gradexact_t;
  Q_T = -std::pow(cos(B_y * y + B_t * t), 0.2e1) * std::pow(cos(D_t * t), 0.2e1) * A_x * A_x * k_1 - std::sin(A_x * x + A_t * t) * std::cos(B_y * y + B_t * t) * std::cos(D_t * t) * A_t * rho * cp_0 - std::pow(cos(A_x * x + A_t * t), 0.2e1) * std::pow(cos(D_t * t), 0.2e1) * B_y * B_y * k_1 - std::cos(A_x * x + A_t * t) * std::sin(B_y * y + B_t * t) * std::cos(D_t * t) * B_t * rho * cp_0 - std::cos(A_x * x + A_t * t) * std::cos(B_y * y + B_t * t) * std::sin(D_t * t) * D_t * rho * cp_0 + 0.3e1 * std::pow(cos(A_x * x + A_t * t), 0.3e1) * std::pow(cos(B_y * y + B_t * t), 0.3e1) * std::pow(cos(D_t * t), 0.3e1) * (A_x * A_x + B_y * B_y) * k_2 + (0.2e1 * A_x * A_x * k_1 - std::sin(A_x * x + A_t * t) * std::cos(B_y * y + B_t * t) * std::cos(D_t * t) * A_t * rho * cp_2 + 0.2e1 * B_y * B_y * k_1 - std::cos(A_x * x + A_t * t) * std::sin(B_y * y + B_t * t) * std::cos(D_t * t) * B_t * rho * cp_2 - std::cos(A_x * x + A_t * t) * std::cos(B_y * y + B_t * t) * std::sin(D_t * t) * D_t * rho * cp_2) * std::pow(cos(A_x * x + A_t * t), 0.2e1) * std::pow(cos(B_y * y + B_t * t), 0.2e1) * std::pow(cos(D_t * t), 0.2e1) + (A_x * A_x * k_0 - 0.2e1 * std::pow(cos(B_y * y + B_t * t), 0.2e1) * std::pow(cos(D_t * t), 0.2e1) * A_x * A_x * k_2 - std::sin(A_x * x + A_t * t) * std::cos(B_y * y + B_t * t) * std::cos(D_t * t) * A_t * rho * cp_1 + B_y * B_y * k_0 - 0.2e1 * std::pow(cos(A_x * x + A_t * t), 0.2e1) * std::pow(cos(D_t * t), 0.2e1) * B_y * B_y * k_2 - std::cos(A_x * x + A_t * t) * std::sin(B_y * y + B_t * t) * std::cos(D_t * t) * B_t * rho * cp_1 - std::cos(A_x * x + A_t * t) * std::cos(B_y * y + B_t * t) * std::sin(D_t * t) * D_t * rho * cp_1) * std::cos(A_x * x + A_t * t) * std::cos(B_y * y + B_t * t) * std::cos(D_t * t);
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
  Scalar rho,
  Scalar k_0,
  Scalar k_1,
  Scalar k_2,
  Scalar cp_0,
  Scalar cp_1,
  Scalar cp_2)
{
  Scalar Q_T;
  //Scalar *gradexact_t;
  Q_T = -std::sin(A_x * x + A_t * t) * std::cos(B_y * y + B_t * t) * std::cos(C_z * z + C_t * t) * std::cos(D_t * t) * rho * cp_0 * A_t - std::cos(A_x * x + A_t * t) * std::sin(B_y * y + B_t * t) * std::cos(C_z * z + C_t * t) * std::cos(D_t * t) * rho * cp_0 * B_t - std::cos(A_x * x + A_t * t) * std::cos(B_y * y + B_t * t) * std::sin(C_z * z + C_t * t) * std::cos(D_t * t) * rho * cp_0 * C_t - std::cos(A_x * x + A_t * t) * std::cos(B_y * y + B_t * t) * std::cos(C_z * z + C_t * t) * std::sin(D_t * t) * rho * cp_0 * D_t - std::pow(cos(B_y * y + B_t * t), 0.2e1) * std::pow(cos(C_z * z + C_t * t), 0.2e1) * std::pow(cos(D_t * t), 0.2e1) * k_1 * A_x * A_x - std::pow(cos(A_x * x + A_t * t), 0.2e1) * std::pow(cos(C_z * z + C_t * t), 0.2e1) * std::pow(cos(D_t * t), 0.2e1) * k_1 * B_y * B_y - std::pow(cos(A_x * x + A_t * t), 0.2e1) * std::pow(cos(B_y * y + B_t * t), 0.2e1) * std::pow(cos(D_t * t), 0.2e1) * k_1 * C_z * C_z + 0.3e1 * std::pow(cos(A_x * x + A_t * t), 0.3e1) * std::pow(cos(B_y * y + B_t * t), 0.3e1) * std::pow(cos(C_z * z + C_t * t), 0.3e1) * std::pow(cos(D_t * t), 0.3e1) * (A_x * A_x + B_y * B_y + C_z * C_z) * k_2 + (-std::sin(A_x * x + A_t * t) * std::cos(B_y * y + B_t * t) * std::cos(C_z * z + C_t * t) * std::cos(D_t * t) * rho * cp_1 * A_t - std::cos(A_x * x + A_t * t) * std::sin(B_y * y + B_t * t) * std::cos(C_z * z + C_t * t) * std::cos(D_t * t) * rho * cp_1 * B_t - std::cos(A_x * x + A_t * t) * std::cos(B_y * y + B_t * t) * std::sin(C_z * z + C_t * t) * std::cos(D_t * t) * rho * cp_1 * C_t - std::cos(A_x * x + A_t * t) * std::cos(B_y * y + B_t * t) * std::cos(C_z * z + C_t * t) * std::sin(D_t * t) * rho * cp_1 * D_t + k_0 * A_x * A_x + k_0 * B_y * B_y + k_0 * C_z * C_z - 0.2e1 * std::pow(cos(B_y * y + B_t * t), 0.2e1) * std::pow(cos(C_z * z + C_t * t), 0.2e1) * std::pow(cos(D_t * t), 0.2e1) * k_2 * A_x * A_x - 0.2e1 * std::pow(cos(A_x * x + A_t * t), 0.2e1) * std::pow(cos(C_z * z + C_t * t), 0.2e1) * std::pow(cos(D_t * t), 0.2e1) * k_2 * B_y * B_y - 0.2e1 * std::pow(cos(A_x * x + A_t * t), 0.2e1) * std::pow(cos(B_y * y + B_t * t), 0.2e1) * std::pow(cos(D_t * t), 0.2e1) * k_2 * C_z * C_z) * std::cos(A_x * x + A_t * t) * std::cos(B_y * y + B_t * t) * std::cos(C_z * z + C_t * t) * std::cos(D_t * t) + (-std::sin(A_x * x + A_t * t) * std::cos(B_y * y + B_t * t) * std::cos(C_z * z + C_t * t) * std::cos(D_t * t) * rho * cp_2 * A_t - std::cos(A_x * x + A_t * t) * std::sin(B_y * y + B_t * t) * std::cos(C_z * z + C_t * t) * std::cos(D_t * t) * rho * cp_2 * B_t - std::cos(A_x * x + A_t * t) * std::cos(B_y * y + B_t * t) * std::sin(C_z * z + C_t * t) * std::cos(D_t * t) * rho * cp_2 * C_t - std::cos(A_x * x + A_t * t) * std::cos(B_y * y + B_t * t) * std::cos(C_z * z + C_t * t) * std::sin(D_t * t) * rho * cp_2 * D_t + 0.2e1 * k_1 * A_x * A_x + 0.2e1 * k_1 * B_y * B_y + 0.2e1 * k_1 * C_z * C_z) * std::pow(cos(A_x * x + A_t * t), 0.2e1) * std::pow(cos(B_y * y + B_t * t), 0.2e1) * std::pow(cos(C_z * z + C_t * t), 0.2e1) * std::pow(cos(D_t * t), 0.2e1);
  return Q_T;
}

int main()
{
  Scalar tfield,tfield2;
  Scalar param=1.2;
  Scalar t= 1;
  Scalar x=.5;
  Scalar y=.4;
  Scalar z=.3;

  /// -----------------------------------------------------------------------
  // initalize 1D
  /// -----------------------------------------------------------------------
  masa_init<Scalar>("temp-test-1d","heateq_1d_unsteady_var");
  
  masa_init_param<Scalar>();

  masa_set_param<Scalar>("A_x",param);
  Scalar A_x=param;

  masa_set_param<Scalar>("A_t",param);
  Scalar A_t=param;

  masa_set_param<Scalar>("D_t",param);
  Scalar D_t=param;

  masa_set_param<Scalar>("cp_0",param);
  Scalar cp_0=param;

  masa_set_param<Scalar>("cp_1",param);
  Scalar cp_1=param;

  masa_set_param<Scalar>("cp_2",param);
  Scalar cp_2=param;

  masa_set_param<Scalar>("rho",param);
  Scalar rho=param;

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
  
  tfield   = masa_eval_source_t<Scalar>(x,t);
  tfield2  = SourceQ_t_1d(x,t,A_x,A_t,D_t,rho,k_0,k_1,k_2,cp_0,cp_1,cp_2);

  tfield=std::abs(tfield-tfield2);

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
  masa_init<Scalar>("temp-test-2d","heateq_2d_unsteady_var");
  masa_init_param<Scalar>();

  masa_set_param<Scalar>("A_x",param);
  A_x=param;

  masa_set_param<Scalar>("A_t",param);
  A_t=param;

  masa_set_param<Scalar>("D_t",param);
  D_t=param;

  masa_set_param<Scalar>("cp_0",param);
  cp_0=param;

  masa_set_param<Scalar>("cp_1",param);
  cp_1=param;

  masa_set_param<Scalar>("cp_2",param);
  cp_2=param;

  masa_set_param<Scalar>("rho",param);
  rho=param;

  masa_set_param<Scalar>("k_0",param);
  k_0=param;

  masa_set_param<Scalar>("k_1",param);
  k_1=param;

  masa_set_param<Scalar>("k_2",param);
  k_2=param;

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
  tfield2   = SourceQ_t_2d(x,y,t,A_x,A_t,B_y,B_t,D_t,rho,k_0,k_1,k_2,cp_0,cp_1,cp_2);

  tfield=std::abs(tfield-tfield2);

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
  masa_init<Scalar>("temp-test-3d","heateq_3d_unsteady_var");
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

  masa_set_param<Scalar>("k_1",param);
  k_1=param;

  masa_set_param<Scalar>("k_2",param);
  k_2=param;

  masa_set_param<Scalar>("cp_0",param);
  cp_0=param;

  masa_set_param<Scalar>("cp_1",param);
  cp_1=param;

  masa_set_param<Scalar>("cp_2",param);
  cp_2=param;

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
  
  tfield   = masa_eval_source_t<Scalar>(x,y,z,t);
  tfield2   = SourceQ_t_3d(x,y,z,t,A_x,A_t,B_y,B_t,C_z,C_t,D_t,rho,k_0,k_1,k_2,cp_0,cp_1,cp_2);

  tfield=std::abs(tfield-tfield2);

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
