// -*-c++-*-
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

#include <masa_internal.h>
using namespace MASA;

template <typename Scalar>
MASA::euler_transient_2d<Scalar>::euler_transient_2d()
{
  this->mmsname = "euler_transient_2d";
  this->dimension = 2;

  this->register_var("rho_0",&rho_0);
  this->register_var("rho_x",&rho_x);
  this->register_var("rho_y",&rho_y);
  this->register_var("rho_t",&rho_t);
  this->register_var("u_0",&u_0);
  this->register_var("u_x",&u_x);
  this->register_var("u_y",&u_y);
  this->register_var("u_t",&u_t);
  this->register_var("v_0",&v_0);
  this->register_var("v_x",&v_x);
  this->register_var("v_y",&v_y);
  this->register_var("v_t",&v_t);
  this->register_var("L",&L);
  this->register_var("p_0",&p_0);
  this->register_var("p_x",&p_x);
  this->register_var("p_y",&p_y);
  this->register_var("p_t",&p_t);
  this->register_var("a_vx",&a_vx);
  this->register_var("a_vy",&a_vy);
  this->register_var("a_vt",&a_vt);
  this->register_var("a_ux",&a_ux);
  this->register_var("a_uy",&a_uy);
  this->register_var("a_ut",&a_ut);
  this->register_var("a_rhox",&a_rhox);
  this->register_var("a_rhoy",&a_rhoy);
  this->register_var("a_rhot",&a_rhot);
  this->register_var("a_py",&a_py);
  this->register_var("a_px",&a_px);
  this->register_var("a_pt",&a_pt);
  this->register_var("Gamma",&Gamma);

  this->init_var();

} // done with constructor

template <typename Scalar>
int MASA::euler_transient_2d<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("rho_0",12);
  err += this->set_var("rho_x",12);
  err += this->set_var("rho_y",12);
  err += this->set_var("rho_t",12);
  err += this->set_var("u_0",12);
  err += this->set_var("u_x",12);
  err += this->set_var("u_y",12);
  err += this->set_var("u_t",12);
  err += this->set_var("v_0",12);
  err += this->set_var("v_x",12);
  err += this->set_var("v_y",12);
  err += this->set_var("v_t",12);
  err += this->set_var("L",12);
  err += this->set_var("p_0",12);
  err += this->set_var("p_x",12);
  err += this->set_var("p_y",12);
  err += this->set_var("p_t",12);
  err += this->set_var("a_vx",12);
  err += this->set_var("a_vy",12);
  err += this->set_var("a_vt",12);
  err += this->set_var("a_ux",12);
  err += this->set_var("a_uy",12);
  err += this->set_var("a_ut",12);
  err += this->set_var("a_rhox",12);
  err += this->set_var("a_rhoy",12);
  err += this->set_var("a_rhot",12);
  err += this->set_var("a_py",12);
  err += this->set_var("a_px",12);
  err += this->set_var("a_pt",12);
  err += this->set_var("Gamma",12);

  return err;

} // done with init_var

// ----------------------------------------
// Source Terms
// ----------------------------------------
#include <math.h>

template <typename Scalar>
Scalar MASA::euler_transient_2d<Scalar>::eval_q_e (Scalar x, Scalar y, Scalar t)
{
  Scalar Q_e_t;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar P;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  P = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_t * cos(a_pt * PI * t / L);
  Q_e_t = -a_px * PI * p_x * Gamma * U * sin(a_px * PI * x / L) / (Gamma - 0.1e1) / L + a_py * PI * p_y * Gamma * V * cos(a_py * PI * y / L) / (Gamma - 0.1e1) / L - a_ut * PI * u_t * RHO * U * sin(a_ut * PI * t / L) / L + a_vt * PI * v_t * RHO * V * cos(a_vt * PI * t / L) / L + (U * U + V * V) * a_rhox * PI * rho_x * U * cos(a_rhox * PI * x / L) / L / 0.2e1 - (U * U + V * V) * a_rhoy * PI * rho_y * V * sin(a_rhoy * PI * y / L) / L / 0.2e1 + (U * U + V * V) * a_rhot * PI * rho_t * cos(a_rhot * PI * t / L) / L / 0.2e1 - a_pt * PI * p_t * sin(a_pt * PI * t / L) / (Gamma - 0.1e1) / L + (0.3e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * U * U / L / 0.2e1 - (a_uy * u_y * sin(a_uy * PI * y / L) + a_vx * v_x * sin(a_vx * PI * x / L)) * PI * RHO * U * V / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.3e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * V * V / L / 0.2e1 + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * Gamma * P / (Gamma - 0.1e1) / L;
  return(Q_e_t);
}

template <typename Scalar>
Scalar MASA::euler_transient_2d<Scalar>::eval_q_u (Scalar x, Scalar y, Scalar t)
{
  Scalar Q_u_t;
  Scalar RHO;
  Scalar U;
  Scalar V;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Q_u_t = a_rhox * PI * rho_x * U * U * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * U * V * sin(a_rhoy * PI * y / L) / L - a_uy * PI * u_y * RHO * V * sin(a_uy * PI * y / L) / L - a_px * PI * p_x * sin(a_px * PI * x / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * U / L + a_rhot * PI * rho_t * U * cos(a_rhot * PI * t / L) / L - a_ut * PI * u_t * RHO * sin(a_ut * PI * t / L) / L;
  return(Q_u_t);
}

template <typename Scalar>
Scalar MASA::euler_transient_2d<Scalar>::eval_q_v (Scalar x, Scalar y, Scalar t)
{
  Scalar Q_v_t;
  Scalar RHO;
  Scalar U;
  Scalar V;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Q_v_t = a_rhox * PI * rho_x * U * V * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * V * sin(a_rhoy * PI * y / L) / L - a_vx * PI * v_x * RHO * U * sin(a_vx * PI * x / L) / L + a_py * PI * p_y * cos(a_py * PI * y / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * V / L + a_rhot * PI * rho_t * V * cos(a_rhot * PI * t / L) / L + a_vt * PI * v_t * RHO * cos(a_vt * PI * t / L) / L;
  return(Q_v_t);
}

template <typename Scalar>
Scalar MASA::euler_transient_2d<Scalar>::eval_q_rho (Scalar x, Scalar y, Scalar t)
{
  Scalar Q_rho_t;
  Scalar RHO;
  Scalar U;
  Scalar V;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Q_rho_t = a_rhox * PI * rho_x * U * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * sin(a_rhoy * PI * y / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO / L + a_rhot * PI * rho_t * cos(a_rhot * PI * t / L) / L;
  return(Q_rho_t);
}


// ----------------------------------------
// Analytical Terms
// ----------------------------------------
template <typename Scalar>
Scalar MASA::euler_transient_2d<Scalar>::eval_exact_rho(Scalar x,Scalar y,Scalar t)
{
  Scalar exact_rho;
  exact_rho = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_t * sin(a_rhot * pi * t / L);
  return exact_rho;
}

template <typename Scalar>
Scalar MASA::euler_transient_2d<Scalar>::eval_exact_p(Scalar x,Scalar y,Scalar t)
{
  Scalar exact_p;
  exact_p = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_t * cos(a_pt * pi * t / L);
  return exact_p;
}

template <typename Scalar>
Scalar MASA::euler_transient_2d<Scalar>::eval_exact_u(Scalar x,Scalar y,Scalar t)
{
  Scalar exact_u;
  exact_u = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_t * cos(a_ut * pi * t / L);
  return exact_u;
}

template <typename Scalar>
Scalar MASA::euler_transient_2d<Scalar>::eval_exact_v(Scalar x,Scalar y,Scalar t)
{
  Scalar exact_v;
  exact_v = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_t * sin(a_vt * pi * t / L);
  return exact_v;
}


// ----------------------------------------
// Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::euler_transient_2d);



//---------------------------------------------------------
// AUTOMASA
// Generated on: 2011-10-28 12:30:36
//---------------------------------------------------------
