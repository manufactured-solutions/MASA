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
// $Id: heat-eq.cpp 18184 2011-03-01 20:09:57Z nick $
//
// radiation.cpp: regression testing radiation solution class
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include<tests.h>
#include <cmath>

using namespace MASA;

typedef double Scalar;

int main()
{

  const Scalar thresh = 5 * std::numeric_limits<Scalar>::epsilon();

  int err = 0;
  Scalar source;
  Scalar exact;
  Scalar x;
  Scalar amp=0;

  //Scalar thresh = 5 * std::numeric_limits<Scalar>::epsilon();
  std::vector<Scalar> vec2;

  // initialize the problem
  err += masa_init<Scalar>("radiation","radiation_integrated_intensity");

  // initialize the default parameters
  err += masa_init_param<Scalar>();

  // intialize the various parameters required for Euler 2D
  // call the sanity check routine 
  // (tests that all variables have been initialized)
  err += masa_sanity_check<Scalar>();

  x = 1;
  source = masa_eval_source_u<Scalar>(x);
  exact  = masa_eval_exact_u<Scalar>(x);

  nancheck(source);
  nancheck(exact);

  // reroute stdout for regressions: TODO remove when logger mechanism
  // is used inside masa.
  freopen("/dev/null","w",stdout);

  exact  = masa_eval_exact_u<Scalar>(-1);
  if(abs(exact+1) > thresh)
    {
      std::cout << "MASA Error: radiation not set properly\n";
      std::cout << abs(exact+1) << std::endl;
      return 1;
    }

  exact  = masa_eval_exact_u<Scalar>(1001.2);
  nancheck(exact);


  masa_get_vec<Scalar>("vec_amp",vec2);
  for(std::vector<Scalar>::iterator it = vec2.begin(); it != vec2.end(); it++)
    {
      amp+= *it;
    }

  //if(std::abs(amp - exact) > thresh)
  //  {
  //    std::cout << "MASA REGRESSION FAILURE\n";
  //    return 1;
  //  }

  return err;

}// end program
