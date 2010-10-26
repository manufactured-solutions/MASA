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

#include <iostream>
#include <stdlib.h>

using namespace MASA;
using namespace std;

template<typename Scalar>
int run_regression()
{

  const Scalar threshold = 1.0e-15; // should be small enough to catch any obvious problems
  const Scalar pi = std::acos(Scalar(-1));

  // solutions
  Scalar ufield,ufield2,ufield3;
  Scalar vfield,vfield2,vfield3;
  Scalar efield,efield2,efield3;
  Scalar rho,rho2,rho3;

  Scalar u_an,u_an2,u_an3;
  Scalar v_an,v_an2,v_an3;
  Scalar p_an,p_an2,p_an3;
  Scalar rho_an,rho_an2,rho_an3;
  
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
      u_an = masa_eval_u_an<Scalar>(x);
      v_an = masa_eval_v_an<Scalar>(x);

      // simple source term check
      ufield = masa_eval_u_source<Scalar>(x);
      vfield = masa_eval_v_source<Scalar>(x);      

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
