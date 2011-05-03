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
// $Author$
// $Id$
//
// axi_euler.cpp: These are the MASA class member functions and constructors
// for the Axisymmetric Euler Equations
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h> 

using namespace MASA;

/* ------------------------------------------------
 *
 *         Axisymmetric EULER EQUATION 1D
 *
 *
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
MASA::axi_euler<Scalar>::axi_euler()
{
  this->mmsname = "axisymmetric_euler";
  this->dimension=2;

  this->register_var("R",&R);
  //this->register_var("k",&k);

  this->register_var("p_0",&p_0);
  this->register_var("p_1",&p_1);
  this->register_var("rho_0",&rho_0);
  this->register_var("rho_1",&rho_1);
  this->register_var("u_1",&u_1);
  this->register_var("w_0",&w_0);
  this->register_var("w_1",&w_1);
  this->register_var("a_pr",&a_pr);
  this->register_var("a_pz",&a_pz);
  this->register_var("a_rhor",&a_rhor);
  this->register_var("a_rhoz",&a_rhoz);
  this->register_var("a_ur",&a_ur);
  this->register_var("a_uz",&a_uz);
  this->register_var("a_wr",&a_wr);
  this->register_var("a_wz",&a_wz);
  this->register_var("L",&L);

  this->register_var("Gamma",&Gamma);
  this->register_var("mu",&mu);

  this->init_var();

}//done with constructor

template <typename Scalar>
int MASA::axi_euler<Scalar>::init_var()
{

  int err = 0;

  // currently randomly generated
  err += this->set_var("R",1.01);

  //err += this->set_var("k",.1791);
  err += this->set_var("p_0",.381);
  err += this->set_var("p_1",.1791);
  err += this->set_var("rho_0",.11);
  err += this->set_var("rho_1",1.9);
  err += this->set_var("u_1",.07);
  err += this->set_var("w_0",.11);
  err += this->set_var("w_1",2.14);
  err += this->set_var("a_pr",1.91);
  err += this->set_var("a_pz",.1711);
  err += this->set_var("a_rhor",1.789);
  err += this->set_var("a_rhoz",4.6);
  err += this->set_var("a_ur",12.01);
  err += this->set_var("a_uz",3.01);
  err += this->set_var("a_wr",1.01);
  err += this->set_var("a_wz",.2);
  err += this->set_var("L",1);

  err += this->set_var("Gamma",.1791);
  err += this->set_var("mu",1.1);

  return err;

} // done with variable initializer


// ----------------------------------------
//   Source Terms
// ----------------------------------------

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_q_rho_u(Scalar r,Scalar z)
{
  Scalar Q_u;
  Q_u = p_1 * std::cos(a_pr * pi * r / L) * std::cos(a_pz * pi * z / L) * a_pr * pi / L - u_1 * u_1 * std::pow(std::sin(a_uz * pi * z / L), Scalar(2)) * std::pow(std::cos(a_ur * pi * r / L) - Scalar(1), Scalar(2)) * rho_1 * std::sin(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L) * a_rhor * pi / L + u_1 * (std::cos(a_ur * pi * r / L) - Scalar(1)) * rho_1 * std::cos(a_rhor * pi * r / L) * std::cos(a_rhoz * pi * z / L) * std::sin(a_uz * pi * z / L) * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * a_rhoz * pi / L - Scalar(2) * u_1 * u_1 * std::pow(std::sin(a_uz * pi * z / L), Scalar(2)) * (std::cos(a_ur * pi * r / L) - Scalar(1)) * std::sin(a_ur * pi * r / L) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * a_ur * pi / L + u_1 * (std::cos(a_ur * pi * r / L) - Scalar(1)) * std::cos(a_uz * pi * z / L) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * a_uz * pi / L + u_1 * (std::cos(a_ur * pi * r / L) - Scalar(1)) * w_1 * std::cos(a_wr * pi * r / L) * std::cos(a_wz * pi * z / L) * std::sin(a_uz * pi * z / L) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * a_wz * pi / L + u_1 * u_1 * std::pow(std::sin(a_uz * pi * z / L), Scalar(2)) * std::pow(std::cos(a_ur * pi * r / L) - Scalar(1), Scalar(2)) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) / r;
  return(Q_u);
}

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_q_rho_w(Scalar r,Scalar z)
{
  Scalar Q_w;
  Q_w = -p_1 * std::sin(a_pr * pi * r / L) * std::sin(a_pz * pi * z / L) * a_pz * pi / L - u_1 * std::sin(a_uz * pi * z / L) * rho_1 * std::sin(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L) * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * (std::cos(a_ur * pi * r / L) - Scalar(1)) * a_rhor * pi / L + std::pow(w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L), Scalar(2)) * rho_1 * std::cos(a_rhor * pi * r / L) * std::cos(a_rhoz * pi * z / L) * a_rhoz * pi / L - u_1 * std::sin(a_uz * pi * z / L) * std::sin(a_ur * pi * r / L) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * a_ur * pi / L - u_1 * std::sin(a_uz * pi * z / L) * w_1 * std::sin(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * (std::cos(a_ur * pi * r / L) - Scalar(1)) * a_wr * pi / L + (Scalar(2) * w_0 + Scalar(2) * w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * w_1 * std::cos(a_wr * pi * r / L) * std::cos(a_wz * pi * z / L) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * a_wz * pi / L + u_1 * std::sin(a_uz * pi * z / L) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * (std::cos(a_ur * pi * r / L) - Scalar(1)) / r;
  return(Q_w);
}

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_q_rho(Scalar r,Scalar z)
{
  Scalar Q_rho;
  Q_rho = -(std::cos(a_ur * pi * r / L) - Scalar(1)) * a_rhor * pi * rho_1 * u_1 * std::sin(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L) * std::sin(a_uz * pi * z / L) / L + (w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L) + w_0) * a_rhoz * pi * rho_1 * std::cos(a_rhor * pi * r / L) * std::cos(a_rhoz * pi * z / L) / L - (rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L) + rho_0) * a_ur * pi * u_1 * std::sin(a_ur * pi * r / L) * std::sin(a_uz * pi * z / L) / L + (rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L) + rho_0) * a_wz * pi * w_1 * std::cos(a_wr * pi * r / L) * std::cos(a_wz * pi * z / L) / L + (rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L) + rho_0) * (std::cos(a_ur * pi * r / L) - Scalar(1)) * u_1 * std::sin(a_uz * pi * z / L) / r;
  return(Q_rho);
}

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_q_rho_e(Scalar r,Scalar z)
{
  Scalar Q_e;
  Q_e = Gamma * std::cos(a_pz * pi * z / L) * std::cos(a_pr * pi * r / L) * p_1 * u_1 * (std::cos(a_ur * pi * r / L) - Scalar(1)) * std::sin(a_uz * pi * z / L) * a_pr * pi / L / (Gamma - Scalar(1)) - Gamma * std::sin(a_pz * pi * z / L) * std::sin(a_pr * pi * r / L) * p_1 * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * a_pz * pi / L / (Gamma - Scalar(1)) - std::sin(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L) * std::sin(a_uz * pi * z / L) * u_1 * (std::cos(a_ur * pi * r / L) - Scalar(1)) * (std::pow(w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L), Scalar(2)) + std::pow(std::sin(a_uz * pi * z / L), Scalar(2)) * std::pow(std::cos(a_ur * pi * r / L) - Scalar(1), Scalar(2)) * u_1 * u_1) * a_rhor * pi * rho_1 / L / Scalar(2) + std::cos(a_rhor * pi * r / L) * std::cos(a_rhoz * pi * z / L) * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * (std::pow(w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L), Scalar(2)) + std::pow(std::sin(a_uz * pi * z / L), Scalar(2)) * std::pow(std::cos(a_ur * pi * r / L) - Scalar(1), Scalar(2)) * u_1 * u_1) * a_rhoz * pi * rho_1 / L / Scalar(2) - std::sin(a_uz * pi * z / L) * std::sin(a_ur * pi * r / L) * (p_0 + p_1 * std::sin(a_pr * pi * r / L) * std::cos(a_pz * pi * z / L)) * a_ur * pi * u_1 * Gamma / L / (Gamma - Scalar(1)) - std::sin(a_uz * pi * z / L) * std::sin(a_ur * pi * r / L) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * (std::pow(w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L), Scalar(2)) + Scalar(3) * std::pow(std::sin(a_uz * pi * z / L), Scalar(2)) * std::pow(std::cos(a_ur * pi * r / L) - Scalar(1), Scalar(2)) * u_1 * u_1) * a_ur * pi * u_1 / L / Scalar(2) + (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * u_1 * u_1 * std::cos(a_uz * pi * z / L) * std::sin(a_uz * pi * z / L) * std::pow(std::cos(a_ur * pi * r / L) - Scalar(1), Scalar(2)) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * a_uz * pi / L - std::sin(a_uz * pi * z / L) * (std::cos(a_ur * pi * r / L) - Scalar(1)) * u_1 * w_1 * std::sin(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * a_wr * pi / L + std::cos(a_wr * pi * r / L) * std::cos(a_wz * pi * z / L) * (p_0 + p_1 * std::sin(a_pr * pi * r / L) * std::cos(a_pz * pi * z / L)) * a_wz * pi * w_1 * Gamma / L / (Gamma - Scalar(1)) + std::cos(a_wr * pi * r / L) * std::cos(a_wz * pi * z / L) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * (Scalar(3) * std::pow(w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L), Scalar(2)) + std::pow(std::sin(a_uz * pi * z / L), Scalar(2)) * std::pow(std::cos(a_ur * pi * r / L) - Scalar(1), Scalar(2)) * u_1 * u_1) * a_wz * pi * w_1 / L / Scalar(2) + std::sin(a_uz * pi * z / L) * (std::cos(a_ur * pi * r / L) - Scalar(1)) * u_1 * (p_0 + p_1 * std::sin(a_pr * pi * r / L) * std::cos(a_pz * pi * z / L)) * Gamma / (Gamma - Scalar(1)) / r + std::sin(a_uz * pi * z / L) * (std::cos(a_ur * pi * r / L) - Scalar(1)) * u_1 * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * (std::pow(w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L), Scalar(2)) + std::pow(std::sin(a_uz * pi * z / L), Scalar(2)) * std::pow(std::cos(a_ur * pi * r / L) - Scalar(1), Scalar(2)) * u_1 * u_1) / r / Scalar(2);
  return(Q_e);
}

// ----------------------------------------
//   Analytical Solutions
// ----------------------------------------

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_exact_u(Scalar r,Scalar z)
{
  Scalar exact_u;
  exact_u = u_1 * (std::cos(a_ur * pi * r / L) - Scalar(1)) * std::sin(a_uz * pi * z / L);
  return exact_u;
}

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_exact_w(Scalar r,Scalar z)
{
  Scalar exact_w;
  exact_w = w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L);
  return exact_w;
}

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_exact_p(Scalar r,Scalar z)
{
  Scalar exact_p;
  exact_p = p_0 + p_1 * std::sin(a_pr * pi * r / L) * std::cos(a_pz * pi * z / L);
  return exact_p;
}

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_exact_rho(Scalar r,Scalar z)
{
  Scalar exact_rho;
  exact_rho = rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L);
  return exact_rho;
}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::axi_euler);
