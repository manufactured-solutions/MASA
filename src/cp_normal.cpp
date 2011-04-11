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
// $Author: 
// $Id: 
//
// cp_normal.cpp: These are the SMASA class member functions and constructors
//                For the conjugate prior normal distribution
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <smasa.h> 

using namespace MASA;

/* ------------------------------------------------
 *
 *  Gaussian Distribution
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
MASA::cp_normal<Scalar>::cp_normal()
{
    this->mmsname = "cp_normal";
    this->dimension=1;

    this->register_var("m",&m);
    this->register_var("sigma",&sigma);
    this->register_var("sigma_d",&sigma_d);

    this->register_var("n",&n);


    this->init_var();

}

template <typename Scalar>
int MASA::cp_normal<Scalar>::init_var()
{

  int err = 0;

  err += this->set_var("m",12);
  err += this->set_var("sigma",12);
  err += this->set_var("sigma_d",12);

  return err;
}

template <typename Scalar>
Scalar MASA::cp_normal<Scalar>::eval_likelyhood(Scalar x)
{
  Scalar likelyhood;
  likelyhood = exp(-(n/(2*pow(sigma_d,2)))*pow((x-x_bar),2));
  return likelyhood;
}

template <typename Scalar>
Scalar MASA::cp_normal<Scalar>::eval_prior(Scalar x)
{
  Scalar prior;
  prior = sqrt(2*pi*pow(sigma,2)) * exp(-(1/(2*pow(sigma,2)))*pow((x-m),2));
  return prior;
}

template <typename Scalar>
Scalar MASA::cp_normal<Scalar>::eval_posterior(Scalar x)
{
  Scalar post;
  Scalar sigmap;
  Scalar mp;
  
  sigmap = sqrt(1/((1/pow(sigma,2)) + (n/pow(sigma_d,2))));
  mp     = pow(sigmap,2) * (m/pow(sigma,2) + (n*x_bar/pow(sigma_d,2)));
  post   = sqrt(2*pi*pow(sigmap,2)) * exp(-(1/(2*pow(sigmap,2)))*pow((x-mp),2));
  return post;
}
