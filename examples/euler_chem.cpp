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
// $Author: nick $
// $Id: euler_example.cpp 17232 2011-02-07 23:35:22Z nick $
//
// euler_example.cpp:
//                    this is an example of the API used for calling 
//                    the 1D euler equation with chemistry
// 
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <examples.h>

using namespace MASA;
using namespace std;

typedef double Scalar;

template<typename Scalar>
Scalar temp_function(Scalar T)
{
  // hackish functional here
  // This is an eyeballed fit (focusing on the 5000K-6000K range) 
  // for the equilibrium constant for N2->N+N dissociation
  Scalar K = exp(4+(T-6000)/500);
  return K;
}

int main()
{
  Scalar efield;
  Scalar ufield;
  Scalar N;
  Scalar Ntwo;

  Scalar exact_u;
  Scalar exact_t;
  Scalar exact_rho;
  Scalar exact_N;
  Scalar exact_Ntwo;

  // error handling
  int err=0;
  
  Scalar x;
  int    nx = 200;  // number of points
  int    lx = 10;     // length
  Scalar dx=double(lx)/double(nx);

  // initialize the problem
  err += masa_init<Scalar>("euler-chemistry-example","euler_chem_1d");

  // call the sanity check routine 
  // (tests that all variables have been initialized)
  err += masa_sanity_check<Scalar>();

  // evaluate source terms over the domain (0<x<1, 0<y<1) 
  for(int i=0;i<nx;i++)
    {
      x=i*dx;

      // evalulate source terms
      ufield = masa_eval_source_rho_u  <Scalar>(x);
      efield = masa_eval_source_rho_e  <Scalar>(x);
      N      = masa_eval_source_rho_N  <Scalar>(x,&temp_function);
      Ntwo   = masa_eval_source_rho_N2 <Scalar>(x,&temp_function);

      // evaluate analytical solution terms
      exact_t    = masa_eval_exact_t     <Scalar>(x);
      exact_u    = masa_eval_exact_u     <Scalar>(x);
      exact_rho  = masa_eval_exact_rho   <Scalar>(x);
      exact_N    = masa_eval_exact_rho_N <Scalar>(x);
      exact_Ntwo = masa_eval_exact_rho_N2<Scalar>(x);

      test(ufield);
      test(efield);
      test(N);
      test(Ntwo);

      test(exact_t);
      test(exact_u);
      test(exact_rho);
      test(exact_N);
      test(exact_Ntwo);
      
    } // done with loop

  return err;

}// done 
