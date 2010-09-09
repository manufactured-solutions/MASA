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
// $Id: euler_example.cpp 12695 2010-08-26 03:47:26Z nick $
//
// euler_example.cpp:
// this is an example of the API used for calling the spelart alamaras 
// (RANS) model
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//


// still under contruction, so will need to fill out more as time goes on
#include <masa.h>
#include <iostream>
#include <fstream>

using namespace MASA;

int main()
{
  // declarations
  double x;
  double tempx;

  double ufield;
  double efield;
  double rho;

  double u_an;
  double v_an;
  double p_an;
  double rho_an;

  //problem size
  double lx,ly;
  double dx,dy;
  int nx,ny;

  // initialize
  nx = 10;  // number of points
  lx=1;     // length

  dx=double(lx/nx);

  // initialize the problem 
  masa_init("spelart-alamaras example","rans_sa");

  // initialize the default parameters
  masa_init_param();

  // intialize the various parameters required for Euler 2D
  // call the sanity check routine 
  // (tests that all variables have been initialized)
  masa_sanity_check();

  // evaluate source terms over the domain (0<x<1, 0<y<1) 
  for(int i=0;i<nx;i++)
      {  
	tempx=i*dx;

	// evaluate source terms
	//masa_eval_u_source  (tempx,tempy,&ufield);
	//masa_eval_e_source  (tempx,tempy,&efield);
	//masa_eval_rho_source(tempx,tempy,&rho);
	
	//evaluate analytical solution
	//masa_eval_u_an        (tempx,tempy,&u_an);
	//masa_eval_p_an        (tempx,tempy,&p_an);
	//masa_eval_rho_an      (tempx,tempy,&rho_an);

      }

}// end program
