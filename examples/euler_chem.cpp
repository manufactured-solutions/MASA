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
// $Id: euler_example.cpp 17232 2011-02-07 23:35:22Z nick $
//
// euler_example.cpp:
// this is an example of the API used for calling the 2D euler equation
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa.h>
#include <iostream>
#include <fstream>

using namespace MASA;

typedef double Scalar;

int main()
{

  int    nx = 10;  // number of points
  int    lx = 1;     // length
  Scalar dx=lx/nx;

  // initialize the problem
  masa_init<Scalar>("euler-chemistry-example","euler_chem_1d");

  // initialize the default parameters
  masa_init_param<Scalar>();

  // intialize the various parameters required for Euler 2D
  // call the sanity check routine 
  // (tests that all variables have been initialized)
  masa_sanity_check<Scalar>();

  // evaluate source terms over the domain (0<x<1, 0<y<1) 
  for(int i=0;i<nx;i++)
    {
      

    } // done with loop

}// done 
