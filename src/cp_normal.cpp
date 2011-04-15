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

    // registering a vector
    this->register_vec("vec_data",vec_data);   
    this->register_var("x_bar",&x_bar);

    this->init_var();

}

template <typename Scalar>
int MASA::cp_normal<Scalar>::init_var()
{
  Scalar av = 0;
  int err = 0;

  err += this->set_var("m",12);
  err += this->set_var("sigma",12);
  err += this->set_var("sigma_d",12);

  // set vector
  vec_data.resize(6);
  for(int it = 0;it<int(vec_data.size());it++)
    {
      vec_data[it]=1;
    }

  for(int it = 0;it<int(vec_data.size());it++)
    {
      av +=vec_data[it];
    }
  av = av / (Scalar)vec_data.size();

  //set x_bar to average of data vector
  err += this->set_var("x_bar",av);
    
  return err;
}

template <typename Scalar>
Scalar MASA::cp_normal<Scalar>::eval_likelyhood(Scalar x)
{
  Scalar likelyhood;
  Scalar av = 0;

  for(int it = 0;it<int(vec_data.size());it++)
    {
      av +=vec_data[it];
    }
  av = av / (Scalar)vec_data.size();
  
  //set x_bar to average of data vector
  this->set_var("x_bar",av);
  
  likelyhood = exp(-(vec_data.size()/(2*pow(sigma_d,2)))*pow((x-x_bar),2));
  return likelyhood;
}


template <typename Scalar>
Scalar MASA::cp_normal<Scalar>::eval_loglikelyhood(Scalar x)
{
  Scalar loglikelyhood;
  Scalar av = 0;

  for(int it = 0;it<int(vec_data.size());it++)
    {
      av +=vec_data[it];
    }
  av = av / (Scalar)vec_data.size();
  
  //set x_bar to average of data vector
  this->set_var("x_bar",av);
  
  loglikelyhood = -(vec_data.size()/(2*pow(sigma_d,2)))*pow((x-x_bar),2);
  return loglikelyhood;
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
  Scalar av = 0;
  Scalar post;
  Scalar sigmap;
  Scalar mp;

  for(int it = 0;it<int(vec_data.size());it++)
    {
      av +=vec_data[it];
    }
  av = av / (Scalar)vec_data.size();

  //set x_bar to average of data vector
  this->set_var("x_bar",av);
  
  sigmap = sqrt(1/((1/pow(sigma,2)) + (Scalar(vec_data.size())/pow(sigma_d,2))));
  mp     = pow(sigmap,2) * (m/pow(sigma,2) + (Scalar(vec_data.size())*x_bar/pow(sigma_d,2)));
  post   = sqrt(2*pi*pow(sigmap,2)) * exp(-(1/(2*pow(sigmap,2)))*pow((x-mp),2));

  return post;
}

template <typename Scalar>
Scalar MASA::cp_normal<Scalar>::factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

template <typename Scalar>
Scalar MASA::cp_normal<Scalar>::eval_cen_mom(int k)
{

  Scalar moment;

  if(k%2 == 0 ) // k is even!
    {
      moment = pow(sigma,k) * (factorial(k) / pow(2,k/2) * factorial(k/2));
    }
  else // k is odd
    {
      moment = 0;
    }

  return moment;

}
// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::cp_normal);
