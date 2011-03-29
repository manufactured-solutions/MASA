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
// $Id: heat.cpp 18162 2011-03-01 05:23:07Z nick $
//
// radiation.cpp: These are the MASA class member functions and constructors
//          For Radiation
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h> 

using namespace MASA;

/* ------------------------------------------------
 *
 *   radiation_integrated_intensity
 *
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
MASA::radiation_integrated_intensity<Scalar>::radiation_integrated_intensity()
{
    this->mmsname = "radiation_integrated_intensity";
    this->dimension=1;

    // registering variables
    this->register_var("no_gauss",&no_gauss);   

    // registering a vector
    this->register_vec("vec_mean",&vec_mean);   
    this->register_vec("vec_amp",&vec_amp);   
    this->register_vec("vec_stdev",&vec_stdev);   
    this->init_var();
  
}//done with constructor

template <typename Scalar>
int MASA::radiation_integrated_intensity<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("no_gauss",5);

  // set size of vectors and set default values
  vec_mean.resize(no_gauss);
  vec_amp.resize(no_gauss);
  vec_stdev.resize(no_gauss);

  err += this->set_vec("vec_mean",&vec_mean);
  err += this->set_vec("vec_amp",&vec_amp);
  err += this->set_vec("vec_stdev",&vec_stdev);

  return err;

}

template <typename Scalar>
Scalar MASA::radiation_integrated_intensity<Scalar>::eval_q_u(Scalar x)
{
  Scalar Q_T;
  Q_T = x;
  return Q_T;
}

template <typename Scalar>
Scalar MASA::radiation_integrated_intensity<Scalar>::eval_exact_u(Scalar x)
{
  Scalar exact_t;
  exact_t = x;
  return exact_t;
}

//int x[5] = {2, 3, 5, 7, 11};
//vector<int> vector1(&x[0], &x[5]);

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::radiation_integrated_intensity);
