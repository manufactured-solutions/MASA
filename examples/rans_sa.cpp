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

typedef double Scalar;

int main()
{
  // declarations
  Scalar x;
  Scalar tempx;

  Scalar ufield;
  Scalar efield;
  Scalar rho;

  Scalar u_an;
  Scalar v_an;
  Scalar p_an;
  Scalar rho_an;

  //problem size
  Scalar lx,ly;
  Scalar dx,dy;
  int nx,ny;

  // initialize
  nx = 10;  // number of points
  lx=1;     // length

  dx=lx/nx;

  // initialize the problem 
  masa_init<Scalar>("spelart-alamaras example","rans_sa");

  // initialize the default parameters
  masa_init_param<Scalar>();

  // intialize the various parameters required for Euler 2D
  // call the sanity check routine 
  // (tests that all variables have been initialized)
  masa_sanity_check<Scalar>();

  // evaluate source terms over the domain (0<x<1, 0<y<1) 
  for(int i=0;i<nx;i++)
      {  
	tempx=i*dx;

	// evaluate source terms
	//masa_eval_u_source<Scalar>  (tempx,tempy,&ufield);
	//masa_eval_e_source<Scalar>  (tempx,tempy,&efield);
	//masa_eval_rho_source<Scalar>(tempx,tempy,&rho);
	
	//evaluate analytical solution
	//masa_eval_u_an<Scalar>        (tempx,tempy,&u_an);
	//masa_eval_p_an<Scalar>        (tempx,tempy,&p_an);
	//masa_eval_rho_an<Scalar>      (tempx,tempy,&rho_an);

      }

}// end program
