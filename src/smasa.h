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
//
// smasa.h: stochastic masa functions
//
// $Id: 
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <masa_internal.h> 
#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <stdlib.h>

// Macro for declaring MASA classes with all supported Scalar types
#define MASA_INSTANTIATE_ALL(my_class) template class my_class<double>; \
                                       template class my_class<long double>

namespace MASA 
{

  // ------------------------------------------------------
  // ---- normal distribution --  normal likelyhood -------
  // ------------------------------------------------------
  template <typename Scalar>
  class cp_normal : public manufactured_solution<Scalar>
  {
    using manufactured_solution<Scalar>::pi;

    Scalar m;    
    Scalar sigma;
    Scalar sigma_d;

    Scalar x_bar; // mean of the data
    std::vector<Scalar> vec_data;

  public:
    cp_normal();
    int init_var();

    Scalar factorial      (int);
    Scalar eval_cen_mom   (int);
    Scalar eval_prior     (Scalar);
    Scalar eval_posterior (Scalar);
    Scalar eval_likelyhood(Scalar);
    Scalar eval_loglikelyhood (Scalar);

  };


}// end masa namespace template
