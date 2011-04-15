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
// $Id: euler2d.cpp 19332 2011-03-31 20:05:16Z nick $
//
// euler2d.cpp :program that tests euler2d from masa against known source term
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <tests.h>

using namespace MASA;
using namespace std;

template<typename Scalar>
int run_regression()
{  
  int err = 0;
  Scalar x = 1;
  Scalar source;
  std::vector<Scalar> data;

  // initialize the problem
  err += masa_init<Scalar>("smasa-example-gaussian","cp_normal");

  // call the sanity check routine 
  // (tests that all variables have been initialized)
  err += masa_sanity_check<Scalar>();

  // grab the default data vector:
  masa_get_vec<Scalar>("vec_data",data);

  // change the data vector to new values:
  masa_set_vec<Scalar>("vec_data",data);

  // evaluate likelyhood, prior, posterior
    for(int j=0;j<nx;j++)
      {  
	tempx=i*dx;
	
	likelyhood   = masa_eval_likelyhood(x);	
	prior        = masa_eval_prior(x);	
	posterior    = masa_eval_posterior(x);
	first_moment = masa_eval_central_moment(1);
	  if(first_moment != 0)
	    {
	      cout << "error in first moment";
	      return 1;
	    }

	test(likelyhood);
	test(prior);
	test(posterior);

      }
}

// queue
int main()
{
  int err=0;

  err += run_regression<double>();
  err += run_regression<long double>();

  return err;
}
