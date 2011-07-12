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
// rans_sa.cpp: this is an example of the API used for 
//              calling the spelart alamaras (RANS) model
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <examples.h>

using namespace MASA;

typedef double Scalar;

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
  err = masa_init<Scalar>("spelart-alamaras example","rans_sa");

  // call the sanity check routine 
  // (tests that all variables have been initialized)
  err = masa_sanity_check<Scalar>();

  // evaluate source terms over the domain (0<x<1)
  for(int i=0;i<nx;i++)
      {  
	tempx=i*dx;

	// evaluate source terms
	ufield = masa_eval_source_u  <Scalar> (tempx);
	vfield = masa_eval_source_v  <Scalar> (tempx);
	
	//evaluate analytical solution
	exact_u = masa_eval_exact_u  <Scalar>  (tempx);
	exact_v = masa_eval_exact_v  <Scalar>  (tempx);

	test(ufield);
	test(vfield);

	test(exact_u);
	test(exact_v);

      }

  return err;

}// end program
