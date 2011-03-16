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
// heat.cpp: These are the MASA class member functions and constructors
//          For Radiation
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h> 

using namespace MASA;

/* ------------------------------------------------
 *
 *   radiation_integrated_intensity
 *
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
MASA::fans_sa_transient_d_finite<Scalar>::fans_sa_transient_d_finite()
{
    this->mmsname = "fans_sa_transient_d_finite";
    this->dimension=1;

  this->register_var("u_0",&u_0);
  this->register_var("u_x",&u_x);
  this->register_var("u_y",&u_y);
  this->register_var("v_0",&v_0);
  this->register_var("v_x",&v_x);
  this->register_var("v_y",&v_y);
  this->register_var("rho_0",&rho_0);
  this->register_var("rho_x",&rho_x);
  this->register_var("rho_y",&rho_y);
  this->register_var("p_0",&p_0);
  this->register_var("p_x",&p_x);
  this->register_var("p_y",&p_y);
  this->register_var("a_px",&a_px);
  this->register_var("a_py",&a_py);
  this->register_var("a_rhox",&a_rhox);
  this->register_var("a_rhoy",&a_rhoy);
  this->register_var("a_ux",&a_ux);
  this->register_var("a_uy",&a_uy);
  this->register_var("a_vx",&a_vx);
  this->register_var("a_vy",&a_vy);
  this->register_var("L",&L);
  this->register_var("Gamma",&Gamma);
  this->register_var("mu",&mu);

  this->register_var("nu_sa_0",&nu_sa_0);
  this->register_var("nu_sa_x",&nu_sa_x);
  this->register_var("nu_sa_y",&nu_sa_y);
  this->register_var("nu_sa_t",&nu_sa_y);

  this->register_var("a_nusax",&a_nusax);
  this->register_var("a_nusay",&a_nusay);
  this->register_var("a_nusat",&a_nusat);

  this->init_var();
  
}//done with constructor

template <typename Scalar>
int MASA::fans_sa_transient_d_finite<Scalar>::init_var()
{
  int err = 0;

  // currently randomly generated
  err += this->set_var("u_0",1.23);
  err += this->set_var("u_x",1.1);
  err += this->set_var("u_y",.08);
  err += this->set_var("v_0",12);
  err += this->set_var("v_x",1.6);
  err += this->set_var("v_y",.67);
  err += this->set_var("rho_0",1.02);
  err += this->set_var("rho_x",7.2);
  err += this->set_var("rho_y",9.8);
  err += this->set_var("p_0",1.2);
  err += this->set_var("p_x",.91);
  err += this->set_var("p_y",.623);
  err += this->set_var("a_px",.165);
  err += this->set_var("a_py",.612);
  err += this->set_var("a_rhox",.627);
  err += this->set_var("a_rhoy",.828);
  err += this->set_var("a_ux",.1987);
  err += this->set_var("a_uy",1.189);
  err += this->set_var("a_vx",1.91);
  err += this->set_var("a_vy",2.901);
  err += this->set_var("Gamma",1.01);
  err += this->set_var("mu",.918);
  err += this->set_var("L",3.02);

  err += this->set_var("nu_sa_0",12.0);
  err += this->set_var("nu_sa_x",12.0);
  err += this->set_var("nu_sa_y",12.0);
  err += this->set_var("nu_sa_t",12.0);

  err += this->set_var("a_nusay",12.0);
  err += this->set_var("a_nusax",12.0);
  err += this->set_var("a_nusat",12.0);

  return err;

}

template <typename Scalar>
Scalar MASA::fans_sa_transient_d_finite<Scalar>::eval_q_u(Scalar x,Scalar y)
{
  Scalar Q_u;
  Q_u = x;
  return Q_u;
}

// ----------------------------------------
//   Manufactured Analytical Solutions
// ----------------------------------------

template <typename Scalar>
Scalar MASA::fans_sa_transient_d_finite<Scalar>::eval_exact_u(Scalar x,Scalar y)
{
  Scalar exact_u;
  exact_u = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  return exact_u;
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_d_finite<Scalar>::eval_exact_v(Scalar x,Scalar y)
{
  Scalar exact_v;
  exact_v = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);
  return exact_v;
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_d_finite<Scalar>::eval_exact_p(Scalar x,Scalar y)
{
  Scalar exact_p;
  exact_p = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L);
  return exact_p;
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_d_finite<Scalar>::eval_exact_rho(Scalar x,Scalar y)
{
  Scalar exact_rho;
  exact_rho = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  return exact_rho;
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_d_finite<Scalar>::eval_exact_nu(Scalar x,Scalar y,Scalar t)
{
  Scalar exact_nu;
  exact_nu = nu_sa_0 + nu_sa_x * cos(a_nusax * pi * x / L) + nu_sa_y * cos(a_nusay * pi * y / L) + nu_sa_t * cos(a_nusat * pi * t / L);
  return exact_nu;
}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::fans_sa_transient_d_finite);
