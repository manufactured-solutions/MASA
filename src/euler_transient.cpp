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
  this->register_var("k",&k);

  this->register_var("u_0",&u_0);
  this->register_var("u_x",&u_x);
  this->register_var("u_t",&u_t);

  this->register_var("rho_0",&rho_0);
  this->register_var("rho_x",&rho_x);
  this->register_var("rho_t",&rho_t);

  this->register_var("p_0",&p_0);
  this->register_var("p_x",&p_x);
  this->register_var("p_t",&p_t);

  this->register_var("a_px",&a_px);
  this->register_var("a_pt",&a_pt);

  this->register_var("a_rhox",&a_rhox);
  this->register_var("a_rhot",&a_rhot);

  this->register_var("a_ux",&a_ux);
  this->register_var("a_ut",&a_ut);

  this->register_var("L",&L);
  this->register_var("Gamma",&Gamma);
  this->register_var("mu",&mu);

  // init defaults
  this->init_var();

}

template <typename Scalar>
int MASA::euler_transient_1d<Scalar>::init_var()
{
  int err = 0;

  // randomly generated
  err += this->set_var("k",1.38);

  err += this->set_var("u_0",14.191);
  err += this->set_var("u_x",1.63);
  err += this->set_var("u_t",1.63);

  err += this->set_var("rho_0",91.5);
  err += this->set_var("rho_x",.13);
  err += this->set_var("rho_t",.13);

  err += this->set_var("p_0",12.1984);
  err += this->set_var("p_x",3.151);
  err += this->set_var("p_t",3.151);

  err += this->set_var("a_px",6.151);
  err += this->set_var("a_pt",6.151);

  err += this->set_var("a_rhox",1.2);
  err += this->set_var("a_rhot",1.2);

  err += this->set_var("a_ux",.03);
  err += this->set_var("a_ut",.03);

  err += this->set_var("L",3.02);
  err += this->set_var("Gamma",16.1);
  err += this->set_var("mu",.091);

  return err;
}


template <typename Scalar>
Scalar MASA::euler_transient_1d<Scalar>::eval_q_rho_u(Scalar x,Scalar t)
{
  Scalar Q_u_t;
  Scalar RHO;
  Scalar U;

  RHO = rho_0 + rho_x * std::sin(a_rhox * pi * x / L) + rho_t * std::sin(a_rhot * pi * t / L);
  U = u_0 + u_x * std::sin(a_ux * pi * x / L) + u_t * std::cos(a_ut * pi * t / L);

  Q_u_t = std::cos(a_rhox * PI * x / L) * a_rhox * PI * rho_x * U * U / L + Scalar(0.2e1) * std::cos(a_ux * PI * x / L) * RHO * a_ux * PI * u_x * U / L + std::cos(a_rhot * PI * t / L) * a_rhot * PI * rho_t * U / L - std::sin(a_ut * PI * t / L) * a_ut * PI * u_t * RHO / L - std::sin(a_px * PI * x / L) * a_px * PI * p_x / L;
  return(Q_u_t);


}

template <typename Scalar>
Scalar MASA::euler_transient_1d<Scalar>::eval_q_rho_e(Scalar x,Scalar t)
{
  Scalar Q_e_t;
  Scalar RHO;
  Scalar U;
  Scalar P;

  RHO = rho_0 + rho_x * std::sin(a_rhox * pi * x / L) + rho_t * std::sin(a_rhot * pi * t / L);
  P = p_0 + p_x * std::cos(a_px * pi * x / L) + p_t * std::cos(a_pt * pi * t / L);
  U = u_0 + u_x * std::sin(a_ux * pi * x / L) + u_t * std::cos(a_ut * pi * t / L);

  Q_e_t = std::cos(a_rhot * PI * t / L) * a_rhot * PI * rho_t * U * U / L / Scalar(0.2e1) - std::sin(a_ut * PI * t / L) * a_ut * PI * u_t * RHO * U / L - std::sin(a_pt * PI * t / L) * a_pt * PI * p_t / (Gamma - Scalar(0.1e1)) / L + std::cos(a_rhox * PI * x / L) * std::pow(U, Scalar(0.3e1)) * a_rhox * PI * rho_x / L / Scalar(0.2e1) + std::cos(a_ux * PI * x / L) * P * a_ux * PI * u_x * Gamma / (Gamma - Scalar(0.1e1)) / L + Scalar(0.3e1) / Scalar(0.2e1) * std::cos(a_ux * PI * x / L) * RHO * U * U * a_ux * PI * u_x / L - std::sin(a_px * PI * x / L) * U * a_px * PI * p_x * Gamma / (Gamma - Scalar(0.1e1)) / L;
  return(Q_e_t);

}

template <typename Scalar>
Scalar MASA::euler_transient_1d<Scalar>::eval_q_rho(Scalar x,Scalar t)
{

  Scalar Q_rho_t;
  Scalar RHO;
  Scalar U;

  RHO = rho_0 + rho_x * std::sin(a_rhox * pi * x / L) + rho_t * std::sin(a_rhot * pi * t / L);
  U = u_0 + u_x * std::sin(a_ux * pi * x / L) + u_t * std::cos(a_ut * pi * t / L);

  Q_rho_t = std::cos(a_rhot * PI * t / L) * a_rhot * PI * rho_t / L + std::cos(a_ux * PI * x / L) * RHO * a_ux * PI * u_x / L + std::cos(a_rhox * PI * x / L) * U * a_rhox * PI * rho_x / L;
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
  exact_u = u_0 + u_x * std::sin(a_ux * pi * x / L) + u_t * std::cos(a_ut * pi * t / L);
  return exact_u;
}

template <typename Scalar>
Scalar MASA::euler_transient_1d<Scalar>::eval_exact_p(Scalar x,Scalar t)
{
  Scalar exact_p;
  exact_p = p_0 + p_x * std::cos(a_px * pi * x / L) + p_t * std::cos(a_pt * pi * t / L);
  return exact_p;
}

template <typename Scalar>
Scalar MASA::euler_transient_1d<Scalar>::eval_exact_rho(Scalar x,Scalar t)
{
  Scalar exact_rho;
  exact_rho = rho_0 + rho_x * std::sin(a_rhox * pi * x / L) + rho_t * std::sin(a_rhot * pi * t / L);
  return exact_rho;
}


// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::euler_transient_1d);
//MASA_INSTANTIATE_ALL(MASA::euler_transient_2d);
//MASA_INSTANTIATE_ALL(MASA::euler_transient_3d);

