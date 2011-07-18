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
// c_laplace_example.cpp: example of the API used for
//                        the 2D laplace's equation
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
#include <examples.h>

int main()
{
  int i,j,k;

  // declarations
  double ffield,phi_field;
  double tempx,tempy;

  //problem size
  double lx,ly;
  double dx,dy;
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
  err += masa_init("laplace example","laplace_2d");

  // evaluate source terms over the domain (0<x<1, 0<y<1) 
  for(i=0;i<nx;i++)
    for(j=0;j<nx;j++)
      {  
	tempx=i*dx;
	tempy=j*dy;

	ffield    = masa_eval_2d_source_f (tempx,tempy);
	phi_field = masa_eval_2d_exact_phi(tempx,tempy);

	test(ffield);
	test(phi_field);
	
      }

  return err;

}// end program

