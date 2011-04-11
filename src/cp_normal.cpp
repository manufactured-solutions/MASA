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


}

template <typename Scalar>
int MASA::cp_normal<Scalar>::init_var()
{

  return 1;
}

template <typename Scalar>
Scalar MASA::cp_normal<Scalar>::eval_likelyhood(Scalar x)
{

  return 1;
}

template <typename Scalar>
Scalar MASA::cp_normal<Scalar>::eval_prior(Scalar x)
{


  return 1;
}

template <typename Scalar>
Scalar MASA::cp_normal<Scalar>::eval_posterior(Scalar x)
{

  return 1;
}
