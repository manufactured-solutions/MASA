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
Scalar MASA::axi_euler<Scalar>::eval_q_u(Scalar r,Scalar z)
{
  Scalar Q_u;
  Q_u = p_1 * cos(a_pr * pi * r / L) * cos(a_pz * pi * z / L) * a_pr * pi / L - u_1 * u_1 * pow(sin(a_uz * pi * z / L), Scalar(2)) * pow(cos(a_ur * pi * r / L) - Scalar(1), Scalar(2)) * rho_1 * sin(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) * a_rhor * pi / L + u_1 * (cos(a_ur * pi * r / L) - Scalar(1)) * rho_1 * cos(a_rhor * pi * r / L) * cos(a_rhoz * pi * z / L) * sin(a_uz * pi * z / L) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * a_rhoz * pi / L - Scalar(2) * u_1 * u_1 * pow(sin(a_uz * pi * z / L), Scalar(2)) * (cos(a_ur * pi * r / L) - Scalar(1)) * sin(a_ur * pi * r / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * a_ur * pi / L + u_1 * (cos(a_ur * pi * r / L) - Scalar(1)) * cos(a_uz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * a_uz * pi / L + u_1 * (cos(a_ur * pi * r / L) - Scalar(1)) * w_1 * cos(a_wr * pi * r / L) * cos(a_wz * pi * z / L) * sin(a_uz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * a_wz * pi / L + u_1 * u_1 * pow(sin(a_uz * pi * z / L), Scalar(2)) * pow(cos(a_ur * pi * r / L) - Scalar(1), Scalar(2)) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) / r;
  return(Q_u);
}

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_q_w(Scalar r,Scalar z)
{
  Scalar Q_w;
  Q_w = -p_1 * sin(a_pr * pi * r / L) * sin(a_pz * pi * z / L) * a_pz * pi / L - u_1 * sin(a_uz * pi * z / L) * rho_1 * sin(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * (cos(a_ur * pi * r / L) - Scalar(1)) * a_rhor * pi / L + pow(w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L), Scalar(2)) * rho_1 * cos(a_rhor * pi * r / L) * cos(a_rhoz * pi * z / L) * a_rhoz * pi / L - u_1 * sin(a_uz * pi * z / L) * sin(a_ur * pi * r / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * a_ur * pi / L - u_1 * sin(a_uz * pi * z / L) * w_1 * sin(a_wr * pi * r / L) * sin(a_wz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (cos(a_ur * pi * r / L) - Scalar(1)) * a_wr * pi / L + (Scalar(2) * w_0 + Scalar(2) * w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * w_1 * cos(a_wr * pi * r / L) * cos(a_wz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * a_wz * pi / L + u_1 * sin(a_uz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * (cos(a_ur * pi * r / L) - Scalar(1)) / r;
  return(Q_w);
}

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_q_rho(Scalar r,Scalar z)
{
  Scalar Q_rho;
  Q_rho = -(cos(a_ur * pi * r / L) - Scalar(1)) * a_rhor * pi * rho_1 * u_1 * sin(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) * sin(a_uz * pi * z / L) / L + (w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L) + w_0) * a_rhoz * pi * rho_1 * cos(a_rhor * pi * r / L) * cos(a_rhoz * pi * z / L) / L - (rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) + rho_0) * a_ur * pi * u_1 * sin(a_ur * pi * r / L) * sin(a_uz * pi * z / L) / L + (rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) + rho_0) * a_wz * pi * w_1 * cos(a_wr * pi * r / L) * cos(a_wz * pi * z / L) / L + (rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) + rho_0) * (cos(a_ur * pi * r / L) - Scalar(1)) * u_1 * sin(a_uz * pi * z / L) / r;
  return(Q_rho);
}

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_q_e(Scalar r,Scalar z)
{
  Scalar Q_e;
  Q_e = Gamma * cos(a_pz * pi * z / L) * cos(a_pr * pi * r / L) * p_1 * u_1 * (cos(a_ur * pi * r / L) - Scalar(1)) * sin(a_uz * pi * z / L) * a_pr * pi / L / (Gamma - Scalar(1)) - Gamma * sin(a_pz * pi * z / L) * sin(a_pr * pi * r / L) * p_1 * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * a_pz * pi / L / (Gamma - Scalar(1)) - sin(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) * sin(a_uz * pi * z / L) * u_1 * (cos(a_ur * pi * r / L) - Scalar(1)) * (pow(w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L), Scalar(2)) + pow(sin(a_uz * pi * z / L), Scalar(2)) * pow(cos(a_ur * pi * r / L) - Scalar(1), Scalar(2)) * u_1 * u_1) * a_rhor * pi * rho_1 / L / Scalar(2) + cos(a_rhor * pi * r / L) * cos(a_rhoz * pi * z / L) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * (pow(w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L), Scalar(2)) + pow(sin(a_uz * pi * z / L), Scalar(2)) * pow(cos(a_ur * pi * r / L) - Scalar(1), Scalar(2)) * u_1 * u_1) * a_rhoz * pi * rho_1 / L / Scalar(2) - sin(a_uz * pi * z / L) * sin(a_ur * pi * r / L) * (p_0 + p_1 * sin(a_pr * pi * r / L) * cos(a_pz * pi * z / L)) * a_ur * pi * u_1 * Gamma / L / (Gamma - Scalar(1)) - sin(a_uz * pi * z / L) * sin(a_ur * pi * r / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (pow(w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L), Scalar(2)) + Scalar(3) * pow(sin(a_uz * pi * z / L), Scalar(2)) * pow(cos(a_ur * pi * r / L) - Scalar(1), Scalar(2)) * u_1 * u_1) * a_ur * pi * u_1 / L / Scalar(2) + (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * u_1 * u_1 * cos(a_uz * pi * z / L) * sin(a_uz * pi * z / L) * pow(cos(a_ur * pi * r / L) - Scalar(1), Scalar(2)) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * a_uz * pi / L - sin(a_uz * pi * z / L) * (cos(a_ur * pi * r / L) - Scalar(1)) * u_1 * w_1 * sin(a_wr * pi * r / L) * sin(a_wz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * a_wr * pi / L + cos(a_wr * pi * r / L) * cos(a_wz * pi * z / L) * (p_0 + p_1 * sin(a_pr * pi * r / L) * cos(a_pz * pi * z / L)) * a_wz * pi * w_1 * Gamma / L / (Gamma - Scalar(1)) + cos(a_wr * pi * r / L) * cos(a_wz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (Scalar(3) * pow(w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L), Scalar(2)) + pow(sin(a_uz * pi * z / L), Scalar(2)) * pow(cos(a_ur * pi * r / L) - Scalar(1), Scalar(2)) * u_1 * u_1) * a_wz * pi * w_1 / L / Scalar(2) + sin(a_uz * pi * z / L) * (cos(a_ur * pi * r / L) - Scalar(1)) * u_1 * (p_0 + p_1 * sin(a_pr * pi * r / L) * cos(a_pz * pi * z / L)) * Gamma / (Gamma - Scalar(1)) / r + sin(a_uz * pi * z / L) * (cos(a_ur * pi * r / L) - Scalar(1)) * u_1 * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (pow(w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L), Scalar(2)) + pow(sin(a_uz * pi * z / L), Scalar(2)) * pow(cos(a_ur * pi * r / L) - Scalar(1), Scalar(2)) * u_1 * u_1) / r / Scalar(2);
  return(Q_e);
}

// ----------------------------------------
//   Analytical Solutions
// ----------------------------------------

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_exact_u(Scalar r,Scalar z)
{
  Scalar u_exact;
  u_exact = u_1 * (cos(a_ur * pi * r / L) - Scalar(1)) * sin(a_uz * pi * z / L);
  return u_exact;
}

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_exact_w(Scalar r,Scalar z)
{
  Scalar w_exact;
  w_exact = w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L);
  return w_exact;
}

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_exact_p(Scalar r,Scalar z)
{
  Scalar p_exact;
  p_exact = p_0 + p_1 * sin(a_pr * pi * r / L) * cos(a_pz * pi * z / L);
  return p_exact;
}

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_exact_rho(Scalar r,Scalar z)
{
  Scalar rho_exact;
  rho_exact = rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L);
  return rho_exact;
}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::axi_euler);
