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
// heat-eq.cpp: example of the API used for the 2D heat equation
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

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
  // declarations
  Scalar tfield;
  Scalar tempx,tempy;

  //problem size
  Scalar lx,ly;
  Scalar dx,dy;
  int nx,ny;

  // error condition
  int err = 0;

  // initialize
  nx = 10;  // number of points
  ny = 10;  
  lx=1;     // length
  ly=1; 

  dx=lx/nx;
  dy=ly/ny;

  // initialize the problem
  err += masa_init<Scalar>("heat equation example","heateq_2d_steady_const");

  // initialize the default parameters
  err += masa_init_param<Scalar>();

  // intialize the various parameters required for Euler 2D
  // call the sanity check routine 
  // (tests that all variables have been initialized)
  err += masa_sanity_check<Scalar>();

  // evaluate source terms over the domain (0<x<1, 0<y<1) 
  for(int i=0;i<nx;i++)
    for(int j=0;j<nx;j++)
      {  
	tempx=i*dx;
	tempy=j*dy;

	tfield = masa_eval_source_t<Scalar>  (tempx,tempy);

	test(tfield);
	
      }

  return err;

}// end program
