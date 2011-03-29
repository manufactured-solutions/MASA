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
// $Author: nick $
// $Id: 
//
// ablation_example.cpp:
// this is an example of the API used for calling the 1d ablation 
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <limits>

using namespace MASA;
using namespace std;

typedef double Scalar;

Scalar MASA_VAR_DEFAULT = -12345.67;
Scalar uninit = -1.33;
Scalar thresh = 5 * numeric_limits<Scalar>::epsilon();

void test(Scalar input)
{
  if(input - MASA_VAR_DEFAULT > thresh)
    {
      exit(1);
    }

  if(input - uninit > thresh)
    {
      exit(1);
    }

}

int main()
{
  // declarations
  Scalar tempx;

  Scalar ufield;
  Scalar efield;
  Scalar cfield;
  Scalar c3field;
  Scalar rho_cfield;
  Scalar rho_c3field;
  Scalar boundary;

  Scalar e_rho_field;
  Scalar rho;

  Scalar exact_u;
  Scalar exact_t;
  Scalar exact_rho;
  Scalar exact_c;
  Scalar exact_c3;

  //problem size
  Scalar lx;
  Scalar dx;
  int nx;

  // initialize
  nx = 10;  // number of points
  lx=1;     // length

  // error handling
  int err=0;

  dx=(Scalar)lx/(Scalar)nx;

  // initialize the problem
  err += masa_init<Scalar>("ablation-example","navierstokes_ablation_1d_steady");

  // test that all variables have been initialized
  err += masa_sanity_check<Scalar>();

  // evaluate source terms over the domain (0<x<1)
  for(int i=0;i<nx;i++)
    {  
      tempx=i*dx;

      /* evaluate source terms */
      ufield = masa_eval_source_rho_u<Scalar>     (tempx);
      e_rho_field = masa_eval_source_rho_e<Scalar>(tempx);
     
      efield  = masa_eval_source_e <Scalar>        (tempx);
      cfield  = masa_eval_source_C <Scalar>        (tempx);
      c3field = masa_eval_source_C3<Scalar>        (tempx);
      rho_cfield  = masa_eval_source_rho_C <Scalar> (tempx);
      rho_c3field = masa_eval_source_rho_C3<Scalar> (tempx);
      boundary   = masa_eval_source_boundary<Scalar> (tempx);
	
      /* evaluate manufactured analytical solution */
      exact_u   = masa_eval_exact_u  <Scalar>   (tempx);
      exact_t   = masa_eval_exact_t  <Scalar>   (tempx);
      exact_rho = masa_eval_exact_rho<Scalar>   (tempx);
      exact_c   = masa_eval_exact_rho_C<Scalar> (tempx);
      exact_c3  = masa_eval_exact_rho_C3<Scalar>(tempx);

      test(ufield);
      test(efield);
      test(cfield);
      test(rho_c3field);
      test(rho_cfield);
      test(rho_c3field);
      test(boundary);
      test(e_rho_field);

      test(exact_u);
      test(exact_t);
      test(exact_rho);
      test(exact_c);
      test(exact_c3);

    }

  return err;

}// end program
