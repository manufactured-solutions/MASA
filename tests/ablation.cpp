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
// $Author:
// $Id:
//
// euler_chem.cpp: program that tests euler with chemistry
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <tests.h>

using namespace MASA;
using namespace std;

// ----------------------------------------
//   Regresssion
// ----------------------------------------

template<typename Scalar>
int run_regression()
{

  Scalar threshold = 5 * numeric_limits<Scalar>::epsilon();

  Scalar ufield,e,erho,c,c3,crho,c3rho;
  Scalar exact_t,exact_u,exact_rho,exact_C,exact_C3;
  Scalar boundary;

  // parameters
  Scalar x;

  //problem size
  int nx = 200;  // number of points
  int lx=10;     // length
  Scalar dx=Scalar(lx)/Scalar(nx);

  // initalize
  masa_init<Scalar>("ablation-test","navierstokes_ablation_1d_steady");

  // evaluate MMS (1D)
  for(int i=0;i<nx;i++)
    {
      x=i*dx;

      // evalulate source terms
      ufield = masa_eval_source_rho_u  <Scalar>(x);
      e      = masa_eval_source_e      <Scalar>(x,&temp_function);
      erho   = masa_eval_source_rho_e  <Scalar>(x);
      c      = masa_eval_source_C      <Scalar>(x);
      c3     = masa_eval_source_C3     <Scalar>(x);
      crho   = masa_eval_source_rho_C  <Scalar>(x);
      c3rho  = masa_eval_source_rho_C3 <Scalar>(x);

      boundary = masa_eval_source_boundary <Scalar>(x);

      // evaluate analytical solution terms
      exact_t    = masa_eval_exact_t     <Scalar>(x);
      exact_u    = masa_eval_exact_u     <Scalar>(x);
      exact_rho  = masa_eval_exact_rho   <Scalar>(x);
      exact_C    = masa_eval_exact_rho_C <Scalar>(x);
      exact_C3   = masa_eval_exact_rho_C3<Scalar>(x);

    } // done w/ spatial interations

  return 0;

} // done with tests

int main()
{
  int err=0;

#ifndef portland_compiler
  err += run_regression<double>();
  err += run_regression<long double>();
#endif

  return err;
}

