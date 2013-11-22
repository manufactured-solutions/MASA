// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012,2013 The PECOS Development Team
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
// c_laplace.c: regression testing laplace solution class
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include<tests.h>

const double thresh = 1.0e-15; // should be small enough to catch any obvious problems

int main()
{

  // variables: 
  double Lx;
  double Ly;
  double field,field2,field3;
  double exact_phi,exact_phi2,exact_phi3;

  // parameters
  double x;
  double y;
  int i,j,k;

  // initalize
  int err;
  int nx = 10;  // number of points
  int ny = 8;  
  int lx=2;     // length
  int ly=1; 
  
  double dx=(double)lx/(double)nx;
  double dy=(double)ly/(double)ny;

  masa_init("laplace regression test","laplace_2d");

  // set params
  masa_init_param();

  // get vars for comparison
  Ly = masa_get_param("Ly");
  Lx = masa_get_param("Lx");

  // check that all terms have been initialized
  err = masa_sanity_check();
  if(err != 0)
    {
      printf( "MASA :: Sanity Check Failed!\n");
      exit(1);
    }  

  // evaluate source terms (2d)
  for(i=0;i<nx;i++)
    for(j=0;j<ny;j++)    
      {
	x=i*dx;
	y=j*dy;

	//evalulate source terms
	field     = masa_eval_2d_source_f  (x,y);

	//evaluate analytical terms
	exact_phi = masa_eval_2d_exact_phi (x,y);

	// get 'exact' solution
	field2  = 2*pow(Lx-x,2) - 8*(Lx-x)*(Lx+x) + 2*pow(Lx+x,2);
	field2 += 2*pow(Ly-y,2) - 8*(Ly-y)*(Ly+y) + 2*pow(Ly+y,2);
	
	exact_phi2 = pow(Ly-y,2)*pow(Ly+y,2) + pow(Lx-x,2)*pow(Lx+x,2);
	
	// test the result is roughly zero
	// choose between abs and rel error
#ifdef MASA_STRICT_REGRESSION

	field3     = fabs(field-field2);
	exact_phi3 = fabs(exact_phi-exact_phi2);

#else

	field3     = fabs(field-field2)/fabs(field2);
	exact_phi3 = fabs(exact_phi-exact_phi2)/fabs(exact_phi2);

#endif	
	threshcheck(field3,thresh);
	threshcheck(exact_phi3,thresh);

      } // done iterating

  // tests passed
  return 0;
  
} // end run_regression
