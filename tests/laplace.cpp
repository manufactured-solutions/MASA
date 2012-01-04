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
// $Author: nick $
// $Id: 
//
// laplace.cpp: regression testing laplace solution class
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include<tests.h>
#include <cmath>

using namespace MASA;

template<typename Scalar>
int run_regression()
{
  const Scalar thresh = 5 * std::numeric_limits<Scalar>::epsilon();

  // variables: 
  Scalar Lx;
  Scalar Ly;
  Scalar field,field2,field3;
  Scalar exact_phi,exact_phi2,exact_phi3;

  // parameters
  Scalar x;
  Scalar y;

  // initalize
  int err;
  int nx = 10;  // number of points
  int ny = 8;  
  int lx=2;     // length
  int ly=1; 
  
  Scalar dx=Scalar(lx)/Scalar(nx);
  Scalar dy=Scalar(ly)/Scalar(ny);

  masa_init<Scalar>("laplace regression test","laplace_2d");

  // set params
  masa_init_param<Scalar>();

  // get vars for comparison
  Ly = masa_get_param<Scalar>("Ly");
  Lx = masa_get_param<Scalar>("Lx");

  // check that all terms have been initialized
  err = masa_sanity_check<Scalar>();
  if(err != 0)
    {
      std::cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }  

  // evaluate source terms (2D)
  for(int i=0;i<nx;i++)
    for(int j=0;j<ny;j++)    
      {
	x=i*dx;
	y=j*dy;

	//evalulate source terms
	field     = masa_eval_source_f  <Scalar>(x,y);

	//evaluate analytical terms
	exact_phi = masa_eval_exact_phi <Scalar>(x,y);

	// get 'exact' solution
	field2  = 2*std::pow(Lx-x,2) - 8*(Lx-x)*(Lx+x) + 2*std::pow(Lx+x,2);
	field2 += 2*std::pow(Ly-y,2) - 8*(Ly-y)*(Ly+y) + 2*std::pow(Ly+y,2);
	
	exact_phi2 = std::pow(Ly-y,2)*std::pow(Ly+y,2) + std::pow(Lx-x,2)*std::pow(Lx+x,2);
	
	// test the result is roughly zero
	// choose between abs and rel error
#ifdef MASA_STRICT_REGRESSION

	field3     = std::abs(field-field2);
	exact_phi3 = std::abs(exact_phi-exact_phi2);

#else

	field3     = std::abs(field-field2)/std::abs(field2);
	exact_phi3 = std::abs(exact_phi-exact_phi2)/std::abs(exact_phi2);

#endif	
	threshcheck(field3,thresh);
	threshcheck(exact_phi3,thresh);

      } // done iterating

  // tests passed
  return 0;
  
} // end run_regression

int main()
{
  int err=0;

  err += run_regression<double>();
  err += run_regression<long double>();

  return err;
}

