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
// $Id: euler_example.cpp 19315 2011-03-31 18:26:04Z nick $
//
// cp_gaussian.cpp: conjugate prior gaussian distribution
// 
// this is an example of the API used for calling the smasa mms
// for a normal prior -- normal likelyhood -- normal posterior
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa.h>

using namespace MASA;

typedef double Scalar;

int main()
{
  int err = 0;
  Scalar x = 1;
  Scalar source,mean,variance;
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
  source = masa_eval_likelyhood(x);
  masa_test_default(source);

  source = masa_eval_prior(x);
  masa_test_default(source);

  source = masa_eval_posterior(x);
  masa_test_default(source);

  // all odd moments must be 0 
  double first_moment  = masa_eval_central_moment<Scalar>(1);
  if(first_moment != 0)
    {
      std::cout << "error in first moment";
      return 1;
    }

  double second_moment = masa_eval_central_moment<Scalar>(2);

  // calculate mean and variance:
  mean = masa_eval_posterior_mean<Scalar>();
  masa_test_default(mean);
  variance = masa_eval_posterior_variance<Scalar>();
  masa_test_default(variance);

  return err;
}
