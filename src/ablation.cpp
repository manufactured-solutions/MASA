// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010 The PECOS Development Team
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
// ablation.cpp: These are the MASA class member functions and constructors
//               For ablation+flow
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h> 

using namespace MASA;

/* ------------------------------------------------
 *   Ablation
 *   Navier Stokes 1D
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
MASA::navierstokes_ablation_1d_steady<Scalar>::navierstokes_ablation_1d_steady()
{
    this->mmsname = "navierstokes_ablation_1d_steady";
    this->dimension=1;

    this->register_var("k",&k);

    this->init_var();
  
}//done with constructor

template <typename Scalar>
int MASA::navierstokes_ablation_1d_steady<Scalar>::init_var()
{
  int err = 0;
  err += this->set_var("k",.82);
  return err;

}

template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_q_rho_u(Scalar x)
{
  Scalar Q_T;
  Q_T = k;
  return Q_T;
}

template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_q_rho_e(Scalar x)
{
  Scalar Q_T;
  Q_T = k;
  return Q_T;  
}

template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_q_rho(Scalar x)
{
  Scalar Q_T;
  Q_T = k;
  return Q_T;  
}

template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_q_rho_C(Scalar x)
{
  Scalar Q_T;
  Q_T = k;
  return Q_T;  
}

template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_q_rho_C3(Scalar x)
{
  Scalar Q_T;
  Q_T = k;
  return Q_T;  
}

template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_q_u_boundary(Scalar x)
{
  Scalar Q_T;
  Q_T = k;
  return Q_T;  
}

/* ------------------------------------------------
 * 
 *
 *    Manufactured Solutions
 * 
 * -----------------------------------------------
 */ 

template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_exact_u(Scalar x)
{
  Scalar Q_T;
  Q_T = k;
  return Q_T;  
}

template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_exact_t(Scalar x)
{
  Scalar Q_T;
  Q_T = k;
  return Q_T;  
}

template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_exact_rho_C(Scalar x)
{
  Scalar Q_T;
  Q_T = k;
  return Q_T;  
}

template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_exact_rho_C3(Scalar x)
{
  Scalar Q_T;
  Q_T = k;
  return Q_T;  
}

MASA_INSTANTIATE_ALL(MASA::navierstokes_ablation_1d_steady);

