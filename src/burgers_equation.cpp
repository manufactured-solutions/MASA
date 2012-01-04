// -*-c++-*-
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012 The PECOS Development Team
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
MASA::burgers_equation<Scalar>::burgers_equation()
{
  this->mmsname = "burgers_equation";
  this->dimension = 2;

  this->register_var("nu",&nu);
  this->register_var("u_0",&u_0);
  this->register_var("u_x",&u_x);
  this->register_var("u_y",&u_y);
  this->register_var("u_t",&u_t);
  this->register_var("v_0",&v_0);
  this->register_var("v_x",&v_x);
  this->register_var("v_y",&v_y);
  this->register_var("v_t",&v_t);
  this->register_var("a_ux",&a_ux);
  this->register_var("a_uy",&a_uy);
  this->register_var("a_ut",&a_ut);
  this->register_var("a_vx",&a_vx);
  this->register_var("a_vy",&a_vy);
  this->register_var("a_vt",&a_vt);
  this->register_var("L",&L);

  this->init_var();

} // done with constructor

template <typename Scalar>
int MASA::burgers_equation<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("nu",12);
  err += this->set_var("u_0",12);
  err += this->set_var("u_x",12);
  err += this->set_var("u_y",12);
  err += this->set_var("u_t",12);
  err += this->set_var("v_0",12);
  err += this->set_var("v_x",12);
  err += this->set_var("v_y",12);
  err += this->set_var("v_t",12);
  err += this->set_var("a_ux",12);
  err += this->set_var("a_uy",12);
  err += this->set_var("a_ut",12);
  err += this->set_var("a_vx",12);
  err += this->set_var("a_vy",12);
  err += this->set_var("a_vt",12);
  err += this->set_var("L",12);

  return err;

} // done with init_var

// ----------------------------------------
// Source Terms
// ----------------------------------------

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_q_v_transient_viscous (Scalar x, Scalar y, Scalar t)
{
  double Qv_tv;
  double U;
  double V;
  double Q_v_time;
  double Q_v_convection;
  double Q_v_dissipation;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / Lt);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / Lt);

  // "Contribution from the time derivative  to the total source term "
  Q_v_time = a_vt * PI * v_t * cos(a_vt * PI * t / Lt) / Lt;

  // "Contribution from the convective terms to the total source term "
  Q_v_convection = -a_vx * PI * v_x * U * sin(a_vx * PI * x / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * V / L;

  // "Contribution from the viscous/dissipation terms to the total source term "
  Q_v_dissipation = a_vx * a_vx * PI * PI * v_x * nu * cos(a_vx * PI * x / L) * pow(L, -0.2e1) + a_vy * a_vy * PI * PI * v_y * nu * sin(a_vy * PI * y / L) * pow(L, -0.2e1);

  // "Total source term "
  Qv_tv = Q_v_dissipation + Q_v_convection + Q_v_time;
  return(Qv_tv);
}

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_q_v_steady_viscous (Scalar x, Scalar y)
{
  double Qv_sv;
  double U;
  double V;
  double Q_v_convection;
  double Q_v_dissipation;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L);

  // "Contribution from the convective terms to the total source term "
  Q_v_convection = -a_vx * PI * v_x * U * sin(a_vx * PI * x / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * V / L;

  // "Contribution from the viscous/dissipation terms to the total source term"
  Q_v_dissipation = a_vx * a_vx * PI * PI * v_x * nu * cos(a_vx * PI * x / L) * pow(L, -0.2e1) + a_vy * a_vy * PI * PI * v_y * nu * sin(a_vy * PI * y / L) * pow(L, -0.2e1);

  // "Total source term "
  Qv_sv = Q_v_dissipation + Q_v_convection;
  return(Qv_sv);
  
}

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_q_v_transient_inviscid (Scalar x, Scalar y, Scalar t)
{
  double Qv_tinv;
  double U;
  double V;
  double Q_v_time;
  double Q_v_convection;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / Lt);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / Lt);

  // "Contribution from the time derivative  to the total source term "  
  Q_v_time = a_vt * PI * v_t * cos(a_vt * PI * t / Lt) / Lt;

  // "Contribution from the convective terms to the total source term "
  Q_v_convection = -a_vx * PI * v_x * U * sin(a_vx * PI * x / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * V / L;

  // "Total source term "
  Qv_tinv = Q_v_convection + Q_v_time;
  return(Qv_tinv);

}

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_q_v_steady_inviscid (Scalar x, Scalar y)
{
  double Qv_sinv;
  double U;
  double V;
  double Q_v_convection;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L);

  // "Contribution from the convective terms to the total source term "
  Q_v_convection = -a_vx * PI * v_x * U * sin(a_vx * PI * x / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * V / L;
  
  // "Total source term "
  Qv_sinv = Q_v_convection;
  return(Qv_sinv);

}

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_q_u_transient_viscous (Scalar x,Scalar y,Scalar t)
{
  double Qu_tv;
  double U;
  double V;
  double Q_u_time;
  double Q_u_convection;
  double Q_u_dissipation;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / Lt);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / Lt);

  // "Contribution from the time derivative  to the total source term"
  Q_u_time = -a_ut * PI * u_t * sin(a_ut * PI * t / Lt) / Lt;

  // "Contribution from the convective terms to the total source term "
  Q_u_convection = -a_uy * PI * u_y * V * sin(a_uy * PI * y / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * U / L;

  // "Contribution from the viscous/dissipation terms to the total source term "
  Q_u_dissipation = a_ux * a_ux * PI * PI * u_x * nu * sin(a_ux * PI * x / L) * pow(L, -0.2e1) + a_uy * a_uy * PI * PI * u_y * nu * cos(a_uy * PI * y / L) * pow(L, -0.2e1);

  // "Total source term"
  Qu_tv = Q_u_dissipation + Q_u_convection + Q_u_time;
  return(Qu_tv);
}

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_q_u_steady_viscous (Scalar x, Scalar y)
{
  double Qu_sv;
  double U;
  double V;
  double Q_u_convection;
  double Q_u_dissipation;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L);

  // "Contribution from the convective terms to the total source term"
  Q_u_convection = -a_uy * PI * u_y * V * sin(a_uy * PI * y / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * U / L;

  // "Contribution from the viscous/dissipation terms to the total source term"
  Q_u_dissipation = a_ux * a_ux * PI * PI * u_x * nu * sin(a_ux * PI * x / L) * pow(L, -0.2e1) + a_uy * a_uy * PI * PI * u_y * nu * cos(a_uy * PI * y / L) * pow(L, -0.2e1);

  // "Total source term"
  Qu_sv = Q_u_dissipation + Q_u_convection;
  return(Qu_sv);
}

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_q_u_transient_inviscid (Scalar x, Scalar y, Scalar t)
{ 
  double Qu_tinv;
  double U;
  double V;
  double Q_u_time;
  double Q_u_convection;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / Lt);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / Lt);

  // "Contribution from the time derivative  to the total source term "
  Q_u_time = -a_ut * PI * u_t * sin(a_ut * PI * t / Lt) / Lt;

  // "Contribution from the convective terms to the total source term "
  Q_u_convection = -a_uy * PI * u_y * V * sin(a_uy * PI * y / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * U / L;

  // "Total source term "
  Qu_tinv = Q_u_convection + Q_u_time;
  return(Qu_tinv);
  
}

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_q_u_steady_inviscid (Scalar x, Scalar y)
{
  double Qu_sinv;
  double U;
  double V;
  double Q_u_convection;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L);

  // "Contribution from the convective terms to the total source term "
  Q_u_convection = -a_uy * PI * u_y * V * sin(a_uy * PI * y / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * U / L;
  
  // "Total source term "
  Qu_sinv = Q_u_convection;
  return(Qu_sinv);
}


// ----------------------------------------
// Analytical Terms
// ----------------------------------------
// public, but will be called from eval_exact_t
template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_exact_u(Scalar x,Scalar y)
{
  Scalar u_an;
  u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  return u_an;
}

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_exact_v(Scalar x,Scalar y)
{
  Scalar v_an;
  v_an = v_0 + v_x * sin(a_vx * pi * x / L) + v_y * cos(a_vy * pi * y / L);
  return v_an;
}

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_exact_u_t(Scalar x,Scalar y, Scalar t)
{
  Scalar u_an;
  u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_t * cos(a_ut * pi * t / L);
  return u_an;
}

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_exact_v_t(Scalar x,Scalar y, Scalar t)
{
  Scalar v_an;
  v_an = v_0 + v_x * sin(a_vx * pi * x / L) + v_y * cos(a_vy * pi * y / L) + v_t * cos(a_vt * pi * t / L);
  return v_an;
}


// ----------------------------------------
// Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::burgers_equation);



//---------------------------------------------------------
// AUTOMASA
// Generated on: 2011-10-26 17:35:23
// updated     : 2011-11-17
//---------------------------------------------------------
