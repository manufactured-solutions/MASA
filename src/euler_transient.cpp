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
// $Id: euler.cpp 17621 2011-02-14 16:53:09Z nick $
//
// euler.cpp: These are the MASA class member functions and constructors
//          For the Euler Equations
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h> 

using namespace MASA;

/* ------------------------------------------------
 *
 *         EULER EQUATION 1D
 *
 *
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
MASA::euler_transient_1d<Scalar>::euler_transient_1d()
{
  this->mmsname = "euler_transient_1d";
  this->dimension=1;

}

template <typename Scalar>
int MASA::euler_transient_1d<Scalar>::init_var()
{
  int err = 0;

  // randomly generated
  err += this->set_var("k",1.38);


  return 0;
}


template <typename Scalar>
Scalar MASA::euler_transient_1d<Scalar>::eval_q_rho_u(Scalar x,Scalar t)
{
  double Q_u_t;
  double RHO;
  double U;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_t * sin(a_rhot * pi * t / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_t * cos(a_ut * pi * t / L);

  Q_u_t = cos(a_rhox * PI * x / L) * a_rhox * PI * rho_x * U * U / L + 0.2e1 * cos(a_ux * PI * x / L) * RHO * a_ux * PI * u_x * U / L + cos(a_rhot * PI * t / L) * a_rhot * PI * rho_t * U / L - sin(a_ut * PI * t / L) * a_ut * PI * u_t * RHO / L - sin(a_px * PI * x / L) * a_px * PI * p_x / L;
  return(Q_u_t);


}

template <typename Scalar>
Scalar MASA::euler_transient_1d<Scalar>::eval_q_rho_e(Scalar x,Scalar t)
{
  double Q_e_t;
  double RHO;
  double U;
  double P;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_t * sin(a_rhot * pi * t / L);
  P = p_0 + p_x * cos(a_px * pi * x / L) + p_t * cos(a_pt * pi * t / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_t * cos(a_ut * pi * t / L);

  Q_e_t = cos(a_rhot * PI * t / L) * a_rhot * PI * rho_t * U * U / L / 0.2e1 - sin(a_ut * PI * t / L) * a_ut * PI * u_t * RHO * U / L - sin(a_pt * PI * t / L) * a_pt * PI * p_t / (Gamma - 0.1e1) / L + cos(a_rhox * PI * x / L) * pow(U, 0.3e1) * a_rhox * PI * rho_x / L / 0.2e1 + cos(a_ux * PI * x / L) * P * a_ux * PI * u_x * Gamma / (Gamma - 0.1e1) / L + 0.3e1 / 0.2e1 * cos(a_ux * PI * x / L) * RHO * U * U * a_ux * PI * u_x / L - sin(a_px * PI * x / L) * U * a_px * PI * p_x * Gamma / (Gamma - 0.1e1) / L;
  return(Q_e_t);

}

template <typename Scalar>
Scalar MASA::euler_transient_1d<Scalar>::eval_q_rho(Scalar x,Scalar t)
{

  double Q_rho_t;
  double RHO;
  double U;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_t * sin(a_rhot * pi * t / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_t * cos(a_ut * pi * t / L);

  Q_rho_t = cos(a_rhot * PI * t / L) * a_rhot * PI * rho_t / L + cos(a_ux * PI * x / L) * RHO * a_ux * PI * u_x / L + cos(a_rhox * PI * x / L) * U * a_rhox * PI * rho_x / L;
  return(Q_rho_t);

}

/* ------------------------------------------------
 *
 *
 *   Analytical terms
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
Scalar MASA::euler_transient_1d<Scalar>::eval_exact_u(Scalar x,Scalar t)
{
  Scalar exact_u;
  exact_u = u_0 + u_x * sin(a_ux * pi * x / L) + u_t * cos(a_ut * pi * t / L);
  return exact_u;
}

template <typename Scalar>
Scalar MASA::euler_transient_1d<Scalar>::eval_exact_p(Scalar x,Scalar t)
{
  Scalar exact_p;
  exact_p = p_0 + p_x * cos(a_px * pi * x / L) + p_t * cos(a_pt * pi * t / L);
  return exact_p;
}

template <typename Scalar>
Scalar MASA::euler_transient_1d<Scalar>::eval_exact_rho(Scalar x,Scalar t)
{
  Scalar exact_rho;
  exact_rho = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_t * sin(a_rhot * pi * t / L);
  return exact_rho;
}


// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::euler_transient_1d);
//MASA_INSTANTIATE_ALL(MASA::euler_transient_2d);
//MASA_INSTANTIATE_ALL(MASA::euler_transient_3d);

