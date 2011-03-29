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
// rans_sa.cpp : program that tests spelart almaras against known source term
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

//
// OBVIOUSLY NOT FINISHED!!
//

#include <config.h>
#include <masa.h>
#include <cmath>
#include <limits>
#include <iostream>
#include <stdlib.h>

using namespace MASA;
using namespace std;

template<typename Scalar>
int run_regression()
{

  // solutions
  Scalar exact_u;
  Scalar exact_v;
  Scalar ufield;
  Scalar vfield;
  
  // parameters
  int nx = 200;  // number of points
  int lx=1;     // length
  Scalar dx=Scalar(lx)/Scalar(nx);

  Scalar x;

  // initalize
  masa_init<Scalar>("spelart alamaras test","rans_sa");

  // initialize the default parameters
  masa_init_param<Scalar>();

  // check that all terms have been initialized
  masa_sanity_check<Scalar>();

  for(int i=0;i<nx;i++)
    {
      x=i*dx;

      // analytical
      exact_u = masa_eval_exact_u<Scalar>(x);
      exact_v = masa_eval_exact_v<Scalar>(x);

      // simple source term check
      ufield = masa_eval_source_u<Scalar>(x);
      vfield = masa_eval_source_v<Scalar>(x);      

    } 
  
  return 0;
}

int main()
{
  int err=0;

  err += run_regression<double>();
  err += run_regression<long double>();

  return err;
}
