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
// $Id: rans_sa.cpp 18184 2011-03-01 20:09:57Z nick $
//
// fans_sa_transient_d_finite:
// this is an example of the API used for calling the spelart alamaras 
// Favre Average Navier Stokes (FANS) model
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

// still under contruction, so will need to fill out more as time goes on
#include <masa.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace MASA;

typedef double Scalar;

Scalar MASA_VAR_DEFAULT = -12345.67;
Scalar uninit = -1.33;

void test(Scalar input)
{
  if(input == MASA_VAR_DEFAULT)
    {
      exit(1);
    }

  if(input == uninit)
    {
      exit(1);
    }

}

int main()
{


  Scalar ufield;
  Scalar vfield;

  // declarations
  Scalar tempx;

  Scalar exact_u;
  Scalar exact_v;

  //error handing
  int err = 0;

  //problem size
  Scalar lx;
  Scalar dx;
  int nx;

  // initialize
  nx = 10;  // number of points
  lx=1;     // length

  dx=lx/nx;

  // initialize the problem 
  err = masa_init<Scalar>("sa example","fans_sa_transient_d_finite");

  // initialize the default parameters
  err = masa_init_param<Scalar>();

  // call the sanity check routine 
  // (tests that all variables have been initialized)
  err = masa_sanity_check<Scalar>();

  return err;

}// end program
