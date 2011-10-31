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
using std::cos;
using std::sin;
using std::pow;
using namespace MASA;

#ifdef portland_compiler
long double pow(long double a, double b){return pow(a,(long double)(b));}
#endif

template <typename Scalar>
MASA::axi_cns_transient<Scalar>::axi_cns_transient()
{
  this->mmsname = "axi_cns_transient";
  this->dimension = 2;

  this->register_var("rho_0",&rho_0);
  this->register_var("rho_r",&rho_r);
  this->register_var("rho_z",&rho_z);
  this->register_var("rho_t",&rho_t);
  this->register_var("a_rhoz",&a_rhoz);
  this->register_var("a_rhor",&a_rhor);
  this->register_var("a_rhot",&a_rhot);
  this->register_var("u_0",&u_0);
  this->register_var("u_r",&u_r);
  this->register_var("u_z",&u_z);
  this->register_var("u_t",&u_t);
  this->register_var("a_uz",&a_uz);
  this->register_var("a_ur",&a_ur);
  this->register_var("a_ut",&a_ut);
  this->register_var("w_0",&w_0);
  this->register_var("w_r",&w_r);
  this->register_var("w_z",&w_z);
  this->register_var("w_t",&w_t);
  this->register_var("a_wz",&a_wz);
  this->register_var("a_wr",&a_wr);
  this->register_var("a_wt",&a_wt);
  this->register_var("p_0",&p_0);
  this->register_var("p_r",&p_r);
  this->register_var("p_z",&p_z);
  this->register_var("p_t",&p_t);
  this->register_var("a_pz",&a_pz);
  this->register_var("a_pr",&a_pr);
  this->register_var("a_pt",&a_pt);
  this->register_var("L",&L);
  this->register_var("Gamma",&Gamma);
  this->register_var("mu",&mu);
  this->register_var("k",&k);
  this->register_var("R",&R);

  this->init_var();

} // done with constructor

template <typename Scalar>
int MASA::axi_cns_transient<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("rho_0",12);
  err += this->set_var("rho_r",12);
  err += this->set_var("rho_z",12);
  err += this->set_var("rho_t",12);
  err += this->set_var("a_rhoz",12);
  err += this->set_var("a_rhor",12);
  err += this->set_var("a_rhot",12);
  err += this->set_var("u_0",12);
  err += this->set_var("u_r",12);
  err += this->set_var("u_z",12);
  err += this->set_var("u_t",12);
  err += this->set_var("a_uz",12);
  err += this->set_var("a_ur",12);
  err += this->set_var("a_ut",12);
  err += this->set_var("w_0",12);
  err += this->set_var("w_r",12);
  err += this->set_var("w_z",12);
  err += this->set_var("w_t",12);
  err += this->set_var("a_wz",12);
  err += this->set_var("a_wr",12);
  err += this->set_var("a_wt",12);
  err += this->set_var("p_0",12);
  err += this->set_var("p_r",12);
  err += this->set_var("p_z",12);
  err += this->set_var("p_t",12);
  err += this->set_var("a_pz",12);
  err += this->set_var("a_pr",12);
  err += this->set_var("a_pt",12);
  err += this->set_var("L",12);
  err += this->set_var("Gamma",12);
  err += this->set_var("mu",12);
  err += this->set_var("k",12);
  err += this->set_var("R",12);

  return err;

} // done with init_var

// ----------------------------------------
// Source Terms
// ----------------------------------------
#include <math.h>

template <typename Scalar>
Scalar MASA::axi_cns_transient<Scalar>::eval_q_e (Scalar r, Scalar z, Scalar t)
{
  Scalar Q_e;
  Scalar RHO;
  Scalar P;
  Scalar U;
  Scalar W;
  RHO = rho_0 + rho_r * cos(a_rhor * PI * r / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_t * sin(a_rhot * PI * t / L);
  P = p_0 + p_r * sin(a_pr * PI * r / L) + p_z * cos(a_pz * PI * z / L) + p_t * cos(a_pt * PI * t / L);
  U = u_r * (cos(a_ur * PI * r / L) - 0.1e1) * (u_z * sin(a_uz * PI * z / L) + u_t * cos(a_ut * PI * t / L));
  W = w_0 + w_r * cos(a_wr * PI * r / L) + w_z * sin(a_wz * PI * z / L) + w_t * cos(a_wt * PI * t / L);
  Q_e = (U * U + W * W) * a_rhot * PI * rho_t * cos(a_rhot * PI * t / L) / L / 0.2e1 - a_pt * PI * p_t * sin(a_pt * PI * t / L) / (Gamma - 0.1e1) / L + (0.3e1 * a_uz * a_uz * u_z * cos(a_ur * PI * r / L) * sin(a_uz * PI * z / L) + 0.4e1 * a_ur * a_ur * u_z * sin(a_uz * PI * z / L) + 0.4e1 * a_ur * a_ur * u_t * cos(a_ut * PI * t / L) - 0.3e1 * a_uz * a_uz * u_z * sin(a_uz * PI * z / L)) * mu * PI * PI * u_r * U * pow(L, -0.2e1) / 0.3e1 - (-a_uz * u_r * u_z * cos(a_ur * PI * r / L) * cos(a_uz * PI * z / L) + a_uz * u_r * u_z * cos(a_uz * PI * z / L) + a_wr * w_r * sin(a_wr * PI * r / L)) * PI * RHO * U * W / L - (a_ur * u_r * u_z * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) + a_ur * u_r * u_t * sin(a_ur * PI * r / L) * cos(a_ut * PI * t / L) - a_wz * w_z * cos(a_wz * PI * z / L)) * Gamma * PI * P / (Gamma - 0.1e1) / L + (sin(a_pr * PI * r / L) * a_pr * a_pr * p_r + cos(a_pz * PI * z / L) * a_pz * a_pz * p_z) * k * PI * PI / R * pow(L, -0.2e1) / RHO - (0.2e1 * a_rhor * a_pr * rho_r * p_r * sin(a_rhor * PI * r / L) * cos(a_pr * PI * r / L) + 0.2e1 * a_rhoz * a_pz * rho_z * p_z * cos(a_rhoz * PI * z / L) * sin(a_pz * PI * z / L)) * k * PI * PI / R * pow(L, -0.2e1) * pow(RHO, -0.2e1) - (U * U + W * W) * a_rhor * PI * rho_r * U * sin(a_rhor * PI * r / L) / L / 0.2e1 + (U * U + W * W) * a_rhoz * PI * rho_z * W * cos(a_rhoz * PI * z / L) / L / 0.2e1 - (a_rhor * a_rhor * rho_r * cos(a_rhor * PI * r / L) + a_rhoz * a_rhoz * rho_z * sin(a_rhoz * PI * z / L)) * k * PI * PI * P / R * pow(L, -0.2e1) * pow(RHO, -0.2e1) - (0.2e1 * a_rhor * a_rhor * rho_r * rho_r * pow(sin(a_rhor * PI * r / L), 0.2e1) + 0.2e1 * a_rhoz * a_rhoz * rho_z * rho_z * pow(cos(a_rhoz * PI * z / L), 0.2e1)) * k * PI * PI * P / R * pow(L, -0.2e1) * pow(RHO, -0.3e1) - a_wt * PI * w_t * RHO * W * sin(a_wt * PI * t / L) / L + RHO * pow(U, 0.3e1) / r / 0.2e1 + 0.4e1 / 0.3e1 * mu * a_ur * a_ur * PI * PI * U * U * pow(L, -0.2e1) - (0.3e1 * a_ur * u_r * u_z * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) + 0.3e1 * a_ur * u_r * u_t * sin(a_ur * PI * r / L) * cos(a_ut * PI * t / L) - a_wz * w_z * cos(a_wz * PI * z / L)) * PI * RHO * U * U / L / 0.2e1 - (a_ur * u_r * u_z * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) + a_ur * u_r * u_t * sin(a_ur * PI * r / L) * cos(a_ut * PI * t / L) - 0.3e1 * a_wz * w_z * cos(a_wz * PI * z / L)) * PI * RHO * W * W / L / 0.2e1 + (a_ur * a_uz * u_r * u_z * sin(a_ur * PI * r / L) * cos(a_uz * PI * z / L) + 0.4e1 * sin(a_wz * PI * z / L) * a_wz * a_wz * w_z) * mu * PI * PI * W * pow(L, -0.2e1) / 0.3e1 + Gamma * P * U / (Gamma - 0.1e1) / r + RHO * U * W * W / r / 0.2e1 - (0.4e1 * a_ur * a_ur * u_r * u_r * u_z * u_z * pow(sin(a_ur * PI * r / L), 0.2e1) * pow(sin(a_uz * PI * z / L), 0.2e1) + 0.8e1 * a_ur * a_ur * u_r * u_r * u_z * u_t * pow(sin(a_ur * PI * r / L), 0.2e1) * sin(a_uz * PI * z / L) * cos(a_ut * PI * t / L) + 0.4e1 * a_ur * a_ur * u_r * u_r * u_t * u_t * pow(sin(a_ur * PI * r / L), 0.2e1) * pow(cos(a_ut * PI * t / L), 0.2e1) + 0.3e1 * a_uz * a_uz * u_r * u_r * u_z * u_z * pow(cos(a_ur * PI * r / L), 0.2e1) * pow(cos(a_uz * PI * z / L), 0.2e1) - 0.6e1 * a_uz * a_uz * u_r * u_r * u_z * u_z * cos(a_ur * PI * r / L) * pow(cos(a_uz * PI * z / L), 0.2e1) + 0.4e1 * a_ur * a_wz * u_r * u_z * w_z * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) * cos(a_wz * PI * z / L) + 0.4e1 * a_ur * a_wz * u_r * u_t * w_z * sin(a_ur * PI * r / L) * cos(a_ut * PI * t / L) * cos(a_wz * PI * z / L) + 0.3e1 * a_uz * a_uz * u_r * u_r * u_z * u_z * pow(cos(a_uz * PI * z / L), 0.2e1) - 0.3e1 * a_uz * a_wr * u_r * u_z * w_r * cos(a_ur * PI * r / L) * cos(a_uz * PI * z / L) * sin(a_wr * PI * r / L) + 0.3e1 * a_uz * a_wr * u_r * u_z * w_r * cos(a_uz * PI * z / L) * sin(a_wr * PI * r / L) + 0.4e1 * a_wz * a_wz * w_z * w_z * pow(cos(a_wz * PI * z / L), 0.2e1)) * mu * PI * PI * pow(L, -0.2e1) / 0.3e1 - Gamma * a_pz * PI * p_z * W * sin(a_pz * PI * z / L) / (Gamma - 0.1e1) / L + 0.4e1 / 0.3e1 * mu * a_wz * PI * w_z * U * cos(a_wz * PI * z / L) / L / r + Gamma * a_pr * PI * p_r * U * cos(a_pr * PI * r / L) / (Gamma - 0.1e1) / L - (cos(a_ur * PI * r / L) - 0.1e1) * a_ut * PI * u_r * u_t * RHO * U * sin(a_ut * PI * t / L) / L - (cos(a_ur * PI * r / L) - 0.1e1) * mu * a_uz * PI * u_r * u_z * W * cos(a_uz * PI * z / L) / L / r / 0.3e1 - k * a_pr * PI * p_r * cos(a_pr * PI * r / L) / R / L / r / RHO - k * a_rhor * PI * rho_r * P * sin(a_rhor * PI * r / L) / R / L / r * pow(RHO, -0.2e1);
  return(Q_e);
}

template <typename Scalar>
Scalar MASA::axi_cns_transient<Scalar>::eval_q_u (Scalar r, Scalar z, Scalar t)
{
  Scalar Q_u;
  Scalar RHO;
  Scalar P;
  Scalar U;
  Scalar W;
  RHO = rho_0 + rho_r * cos(a_rhor * PI * r / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_t * sin(a_rhot * PI * t / L);
  P = p_0 + p_r * sin(a_pr * PI * r / L) + p_z * cos(a_pz * PI * z / L) + p_t * cos(a_pt * PI * t / L);
  U = u_r * (cos(a_ur * PI * r / L) - 0.1e1) * (u_z * sin(a_uz * PI * z / L) + u_t * cos(a_ut * PI * t / L));
  W = w_0 + w_r * cos(a_wr * PI * r / L) + w_z * sin(a_wz * PI * z / L) + w_t * cos(a_wt * PI * t / L);
  Q_u = (cos(a_ur * PI * r / L) - 0.1e1) * a_uz * PI * u_r * u_z * RHO * W * cos(a_uz * PI * z / L) / L - a_rhor * PI * rho_r * U * U * sin(a_rhor * PI * r / L) / L + a_rhoz * PI * rho_z * U * W * cos(a_rhoz * PI * z / L) / L - (cos(a_ur * PI * r / L) - 0.1e1) * a_ut * PI * u_r * u_t * RHO * sin(a_ut * PI * t / L) / L + a_rhot * PI * rho_t * U * cos(a_rhot * PI * t / L) / L + 0.4e1 / 0.3e1 * mu * a_ur * a_ur * PI * PI * U * pow(L, -0.2e1) + a_pr * PI * p_r * cos(a_pr * PI * r / L) / L + (0.3e1 * a_uz * a_uz * u_z * cos(a_ur * PI * r / L) * sin(a_uz * PI * z / L) + 0.4e1 * a_ur * a_ur * u_z * sin(a_uz * PI * z / L) + 0.4e1 * a_ur * a_ur * u_t * cos(a_ut * PI * t / L) - 0.3e1 * a_uz * a_uz * u_z * sin(a_uz * PI * z / L)) * mu * PI * PI * u_r * pow(L, -0.2e1) / 0.3e1 - (0.2e1 * a_ur * u_r * u_z * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) + 0.2e1 * a_ur * u_r * u_t * sin(a_ur * PI * r / L) * cos(a_ut * PI * t / L) - a_wz * w_z * cos(a_wz * PI * z / L)) * PI * RHO * U / L + RHO * U * U / r + 0.2e1 / 0.3e1 * (a_ur * u_r * u_z * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) + a_ur * u_r * u_t * sin(a_ur * PI * r / L) * cos(a_ut * PI * t / L) + a_wz * w_z * cos(a_wz * PI * z / L)) * mu * PI / L / r;
  return(Q_u);
}

template <typename Scalar>
Scalar MASA::axi_cns_transient<Scalar>::eval_q_w (Scalar r, Scalar z, Scalar t)
{
  Scalar Q_w;
  Scalar RHO;
  Scalar P;
  Scalar U;
  Scalar W;
  RHO = rho_0 + rho_r * cos(a_rhor * PI * r / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_t * sin(a_rhot * PI * t / L);
  P = p_0 + p_r * sin(a_pr * PI * r / L) + p_z * cos(a_pz * PI * z / L) + p_t * cos(a_pt * PI * t / L);
  U = u_r * (cos(a_ur * PI * r / L) - 0.1e1) * (u_z * sin(a_uz * PI * z / L) + u_t * cos(a_ut * PI * t / L));
  W = w_0 + w_r * cos(a_wr * PI * r / L) + w_z * sin(a_wz * PI * z / L) + w_t * cos(a_wt * PI * t / L);
  Q_w = -a_rhor * PI * rho_r * U * W * sin(a_rhor * PI * r / L) / L + a_rhoz * PI * rho_z * W * W * cos(a_rhoz * PI * z / L) / L - a_wr * PI * w_r * RHO * U * sin(a_wr * PI * r / L) / L + a_rhot * PI * rho_t * W * cos(a_rhot * PI * t / L) / L - a_wt * PI * w_t * RHO * sin(a_wt * PI * t / L) / L - a_pz * PI * p_z * sin(a_pz * PI * z / L) / L - (cos(a_ur * PI * r / L) - 0.1e1) * mu * a_uz * PI * u_r * u_z * cos(a_uz * PI * z / L) / L / r / 0.3e1 - (a_ur * u_r * u_z * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) + a_ur * u_r * u_t * sin(a_ur * PI * r / L) * cos(a_ut * PI * t / L) - 0.2e1 * a_wz * w_z * cos(a_wz * PI * z / L)) * PI * RHO * W / L + (a_ur * a_uz * u_r * u_z * sin(a_ur * PI * r / L) * cos(a_uz * PI * z / L) + 0.4e1 * a_wz * a_wz * w_z * sin(a_wz * PI * z / L)) * mu * PI * PI * pow(L, -0.2e1) / 0.3e1 + RHO * U * W / r;
  return(Q_w);
}

template <typename Scalar>
Scalar MASA::axi_cns_transient<Scalar>::eval_q_rho (Scalar r, Scalar z, Scalar t)
{
  Scalar Q_rho;
  Scalar RHO;
  Scalar P;
  Scalar U;
  Scalar W;
  RHO = rho_0 + rho_r * cos(a_rhor * PI * r / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_t * sin(a_rhot * PI * t / L);
  P = p_0 + p_r * sin(a_pr * PI * r / L) + p_z * cos(a_pz * PI * z / L) + p_t * cos(a_pt * PI * t / L);
  U = u_r * (cos(a_ur * PI * r / L) - 0.1e1) * (u_z * sin(a_uz * PI * z / L) + u_t * cos(a_ut * PI * t / L));
  W = w_0 + w_r * cos(a_wr * PI * r / L) + w_z * sin(a_wz * PI * z / L) + w_t * cos(a_wt * PI * t / L);
  Q_rho = -a_rhor * PI * rho_r * U * sin(a_rhor * PI * r / L) / L + a_rhoz * PI * rho_z * W * cos(a_rhoz * PI * z / L) / L + a_rhot * PI * rho_t * cos(a_rhot * PI * t / L) / L - (a_ur * u_r * u_z * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) + a_ur * u_r * u_t * sin(a_ur * PI * r / L) * cos(a_ut * PI * t / L) - a_wz * w_z * cos(a_wz * PI * z / L)) * PI * RHO / L + RHO * U / r;
  return(Q_rho);
}



// ----------------------------------------
// Analytical Terms
// ----------------------------------------

template <typename Scalar>
Scalar MASA::axi_cns_transient<Scalar>::eval_exact_rho(Scalar r, Scalar z, Scalar t)
{
  Scalar rho_an;
  rho_an = rho_0 + rho_r * cos(a_rhor * PI * r / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_t * sin(a_rhot * PI * t / L);
  return rho_an;
}

template <typename Scalar>
Scalar MASA::axi_cns_transient<Scalar>::eval_exact_p(Scalar r, Scalar z, Scalar t)
{
  Scalar p_an;
  p_an = p_0 + p_r * sin(a_pr * PI * r / L) + p_z * cos(a_pz * PI * z / L) + p_t * cos(a_pt * PI * t / L);
  return p_an;
}

template <typename Scalar>
Scalar MASA::axi_cns_transient<Scalar>::eval_exact_u(Scalar r, Scalar z, Scalar t)
{
  Scalar u_an;
  u_an = u_r * (cos(a_ur * PI * r / L) - 0.1e1) * (u_z * sin(a_uz * PI * z / L) + u_t * cos(a_ut * PI * t / L));
  return u_an;
}

template <typename Scalar>
Scalar MASA::axi_cns_transient<Scalar>::eval_exact_w(Scalar r, Scalar z, Scalar t)
{
  Scalar w_an;
  w_an = w_0 + w_r * cos(a_wr * PI * r / L) + w_z * sin(a_wz * PI * z / L) + w_t * cos(a_wt * PI * t / L);
  return w_an;
}




// ----------------------------------------
// Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::axi_cns_transient);



//---------------------------------------------------------
// AUTOMASA
// Generated on: 2011-10-31 13:07:51
//---------------------------------------------------------
