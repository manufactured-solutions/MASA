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
// $Id: 
//
// laplace_example.cpp: example of the API used for
//                      the 2D laplace's equation
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
#include <examples.h>

using namespace MASA;
typedef double Scalar;

int main()
{
  // declarations
  Scalar ffield,phi_field;
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
  err += masa_init<Scalar>("laplace example","laplace_2d");

  // call the sanity check routine 
  // (tests that all variables have been initialized)
  err += masa_sanity_check<Scalar>();

  // evaluate source terms over the domain (0<x<1, 0<y<1) 
  for(int i=0;i<nx;i++)
    for(int j=0;j<nx;j++)
      {  
	tempx=i*dx;
	tempy=j*dy;

	ffield    = masa_eval_source_f <Scalar>  (tempx,tempy);
	phi_field = masa_eval_exact_phi<Scalar>  (tempx,tempy);

	test(ffield);
	test(phi_field);
	
      }

  return err;

}// end program

