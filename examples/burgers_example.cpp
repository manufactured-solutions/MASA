// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012 The PECOS Development Team
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
// burgers_example.cpp:
// this is an example of the API used for calling the 2D burgers equation
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa.h>

using namespace MASA;

typedef double Scalar;

int main()
{
  // declarations
  Scalar tempx,tempy, tempt;

  Scalar ufield;
  Scalar vfield;
  Scalar efield;
  Scalar rho;

  Scalar exact_u;
  Scalar exact_v;
  Scalar exact_p;
  Scalar exact_rho;

  //problem size
  Scalar lx,ly;
  Scalar dx,dy;
  int nx,ny;

  // initialize
  nx = 10;  // number of points
  ny = 10;  
  lx=1;     // length
  ly=1; 

  // error handling
  int err=0;

  dx=(Scalar)lx/(Scalar)nx;
  dy=(Scalar)ly/(Scalar)ny;

  // initialize the problem
  err += masa_init<Scalar>("burgers-example","burgers_equation");

  // call the sanity check routine 
  // (tests that all variables have been initialized)
  err += masa_sanity_check<Scalar>();

  tempt = 0.0;

  // evaluate source terms over the domain (0<x<1, 0<y<1) 
  for(int i=0;i<nx;i++)
    for(int j=0;j<nx;j++)
      {  
	tempx=i*dx;
	tempy=j*dy;

	// evaluate source terms
	ufield = masa_eval_source_u<Scalar>  (tempx,tempy,tempt);
	vfield = masa_eval_source_v<Scalar>  (tempx,tempy,tempt);
	
	//evaluate analytical solution
	exact_u   = masa_eval_exact_u  <Scalar>   (tempx,tempy);
	exact_v   = masa_eval_exact_v  <Scalar>   (tempx,tempy);

	masa_test_default(ufield);
	masa_test_default(vfield);

	masa_test_default(exact_u);
	masa_test_default(exact_v);
      }

  return err;

}// end program
