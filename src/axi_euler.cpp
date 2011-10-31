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

using std::cos;
using std::sin;
using std::pow;
using namespace MASA;

// overload for pgi compilers
#ifdef portland_compiler
long double pow(long double a, double b){return pow(a,(long double)(b));}
#endif

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

  this->register_var("p_0",&p_0);
  this->register_var("p_r",&p_r);
  this->register_var("p_z",&p_z);
  this->register_var("rho_0",&rho_0);
  this->register_var("rho_r",&rho_r);
  this->register_var("rho_z",&rho_z);
  this->register_var("u_r",&u_r);
  this->register_var("u_z",&u_z);
  this->register_var("w_0",&w_0);
  this->register_var("w_r",&w_r);
  this->register_var("w_z",&w_z);
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
  err += this->set_var("p_r",.1791);
  err += this->set_var("p_z",.1791);
  err += this->set_var("rho_0",.11);
  err += this->set_var("rho_r",1.9);
  err += this->set_var("rho_z",1.9);
  err += this->set_var("u_r",.07);
  err += this->set_var("u_z",.07);
  err += this->set_var("w_0",.11);
  err += this->set_var("w_r",.11);
  err += this->set_var("w_z",2.14);
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
  Scalar RHO;
  Scalar U;
  Scalar W;
  RHO = rho_0 + rho_r * cos(a_rhor * PI * r / L) + rho_z * sin(a_rhoz * PI * z / L);
  U = u_r * u_z * (cos(a_ur * PI * r / L) - 0.1e1) * sin(a_uz * PI * z / L);
  W = w_0 + w_r * cos(a_wr * PI * r / L) + w_z * sin(a_wz * PI * z / L);
  Q_u = (cos(a_ur * PI * r / L) - 0.1e1) * a_uz * PI * u_r * u_z * RHO * W * cos(a_uz * PI * z / L) / L - a_rhor * PI * rho_r * U * U * sin(a_rhor * PI * r / L) / L + a_rhoz * PI * rho_z * U * W * cos(a_rhoz * PI * z / L) / L + a_pr * PI * p_r * cos(a_pr * PI * r / L) / L - (0.2e1 * a_ur * u_r * u_z * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) - a_wz * w_z * cos(a_wz * PI * z / L)) * PI * RHO * U / L + RHO * U * U / r;
  return(Q_u);
}

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_q_rho_w(Scalar r,Scalar z)
{
  Scalar Q_w;
  Scalar RHO;
  Scalar U;
  Scalar W;
  RHO = rho_0 + rho_r * cos(a_rhor * PI * r / L) + rho_z * sin(a_rhoz * PI * z / L);
  U = u_r * u_z * (cos(a_ur * PI * r / L) - 0.1e1) * sin(a_uz * PI * z / L);
  W = w_0 + w_r * cos(a_wr * PI * r / L) + w_z * sin(a_wz * PI * z / L);
  Q_w = -a_rhor * PI * rho_r * U * W * sin(a_rhor * PI * r / L) / L + a_rhoz * PI * rho_z * W * W * cos(a_rhoz * PI * z / L) / L - a_wr * PI * w_r * RHO * U * sin(a_wr * PI * r / L) / L - a_pz * PI * p_z * sin(a_pz * PI * z / L) / L - (a_ur * u_r * u_z * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) - 0.2e1 * a_wz * w_z * cos(a_wz * PI * z / L)) * PI * RHO * W / L + RHO * U * W / r;
  return(Q_w);
}

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_q_rho(Scalar r,Scalar z)
{
  Scalar Q_rho;
  Scalar RHO;
  Scalar U;
  Scalar W;
  RHO = rho_0 + rho_r * cos(a_rhor * PI * r / L) + rho_z * sin(a_rhoz * PI * z / L);
  U = u_r * u_z * (cos(a_ur * PI * r / L) - 0.1e1) * sin(a_uz * PI * z / L);
  W = w_0 + w_r * cos(a_wr * PI * r / L) + w_z * sin(a_wz * PI * z / L);
  Q_rho = -a_rhor * PI * rho_r * U * sin(a_rhor * PI * r / L) / L + a_rhoz * PI * rho_z * W * cos(a_rhoz * PI * z / L) / L - (a_ur * u_r * u_z * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) - a_wz * w_z * cos(a_wz * PI * z / L)) * PI * RHO / L + RHO * U / r;
  return(Q_rho);
}

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_q_rho_e(Scalar r,Scalar z)
{
  Scalar Q_e;
  Scalar RHO;
  Scalar P;
  Scalar U;
  Scalar W;
  RHO = rho_0 + rho_r * cos(a_rhor * PI * r / L) + rho_z * sin(a_rhoz * PI * z / L);
  P = p_0 + p_r * sin(a_pr * PI * r / L) + p_z * cos(a_pz * PI * z / L);
  U = u_r * u_z * (cos(a_ur * PI * r / L) - 0.1e1) * sin(a_uz * PI * z / L);
  W = w_0 + w_r * cos(a_wr * PI * r / L) + w_z * sin(a_wz * PI * z / L);
  Q_e = -Gamma * a_ur * PI * u_r * u_z * P * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) / (Gamma - 0.1e1) / L + (cos(a_ur * PI * r / L) - 0.1e1) * a_uz * PI * u_r * u_z * RHO * U * W * cos(a_uz * PI * z / L) / L - (0.3e1 * U * U + W * W) * a_ur * PI * u_r * u_z * RHO * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) / L / 0.2e1 - a_wr * PI * w_r * RHO * U * W * sin(a_wr * PI * r / L) / L + Gamma * a_pr * PI * p_r * U * cos(a_pr * PI * r / L) / (Gamma - 0.1e1) / L - Gamma * a_pz * PI * p_z * W * sin(a_pz * PI * z / L) / (Gamma - 0.1e1) / L + Gamma * a_wz * PI * w_z * P * cos(a_wz * PI * z / L) / (Gamma - 0.1e1) / L - (U * U + W * W) * a_rhor * PI * rho_r * U * sin(a_rhor * PI * r / L) / L / 0.2e1 + (U * U + W * W) * a_rhoz * PI * rho_z * W * cos(a_rhoz * PI * z / L) / L / 0.2e1 + (U * U + 0.3e1 * W * W) * a_wz * PI * w_z * RHO * cos(a_wz * PI * z / L) / L / 0.2e1 + Gamma * P * U / (Gamma - 0.1e1) / r + (U * U + W * W) * RHO * U / r / 0.2e1;
  return(Q_e);
}

// ----------------------------------------
//   Analytical Solutions
// ----------------------------------------

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_exact_u(Scalar r,Scalar z)
{
  Scalar u_an;
  u_an = u_r * u_z * (cos(a_ur * PI * r / L) - 0.1e1) * sin(a_uz * PI * z / L);
  return u_an;
}

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_exact_w(Scalar r,Scalar z)
{
  Scalar w_an;
  w_an = w_0 + w_r * cos(a_wr * PI * r / L) + w_z * sin(a_wz * PI * z / L);
  return w_an;
}

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_exact_p(Scalar r,Scalar z)
{
  Scalar p_an;
  p_an = p_0 + p_r * sin(a_pr * PI * r / L) + p_z * cos(a_pz * PI * z / L);
  return p_an;
}

template <typename Scalar>
Scalar MASA::axi_euler<Scalar>::eval_exact_rho(Scalar r,Scalar z)
{
  Scalar rho_an;
  rho_an = rho_0 + rho_r * cos(a_rhor * PI * r / L) + rho_z * sin(a_rhoz * PI * z / L);
  return rho_an;
}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::axi_euler);
