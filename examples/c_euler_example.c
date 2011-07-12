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
// c_euler_example.c: this is an example of the euler equation in 1d for C
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <examples.h>

int main()
{
  double ufield;
  double efield;
  double rhofield;

  double exact_u;
  double exact_p;
  double exact_rho;

  //problem size
  double lx;
  double dx;
  int nx;

  // initialize
  nx = 10;  // number of points
  lx=1;     // length

  // error handling
  int err=0;
  int i;
  double tempx;

  dx=(double)lx/(double)nx;
  
  // init
  err += masa_init("nick","euler_1d");  

  //check all initialized properly
  err += masa_sanity_check();

  for(i=0;i<nx;i++)
    {
      tempx=i*dx;

      ufield   = masa_eval_1d_source_rho_u(tempx);
      efield   = masa_eval_1d_source_rho_e(tempx);
      rhofield = masa_eval_1d_source_rho(tempx);

      //evaluate analytical terms
      exact_u   = masa_eval_1d_exact_u      (tempx);
      exact_p   = masa_eval_1d_exact_p      (tempx);
      exact_rho = masa_eval_1d_exact_rho    (tempx);

      test(ufield);
      test(efield);
      test(rhofield);

      test(exact_u);
      test(exact_p);
      test(exact_rho);

    }

  return err; // done
}
