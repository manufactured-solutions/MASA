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
// $Id: euler_example.cpp 17232 2011-02-07 23:35:22Z nick $
//
// euler_transient.cpp: this is an example of the API used for 
//                      calling the 1d transient euler equation
// 
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa.h>

using namespace MASA;
typedef double Scalar;

int main()
{

  Scalar ufield;
  Scalar efield;
  Scalar rhofield;

  // parameters
  Scalar x;
  Scalar t;

  // error handling
  int err = 0;

  //problem size
  int nx = 200;  // number of points
  int lx=10;     // length
  Scalar dx=Scalar(lx)/Scalar(nx);

  int nt = 100;  // number of points
  int lt = 10;     // length 
  Scalar dt = Scalar(lt)/Scalar(nt);

  // initalize
  err += masa_init<Scalar>("euler-test","euler_transient_1d");

  // check that all terms have been initialized
  err += masa_sanity_check<Scalar>();

  // evaluate MMS (1D)
  for(int i=0;i<nx;i++)
    for(int j=0;j<nt;j++)
      {

	x=i*dx;
	t=j*dt;
	
	// evalulate source terms
	ufield = masa_eval_source_rho_u  <Scalar>(x,t);
	efield = masa_eval_source_rho_e  <Scalar>(x,t);
	rhofield = masa_eval_source_rho  <Scalar>(x,t);

	masa_test_default(ufield);
	masa_test_default(efield);
	masa_test_default(rhofield);
	
      }
  
  return err;

}
