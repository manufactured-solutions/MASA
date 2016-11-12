// -*-c++-*-
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012,2013 The PECOS Development Team
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

#ifdef HAVE_METAPHYSICL

#include <ad_masa.h>

typedef ShadowNumber<double, long double> RawScalar;
const unsigned int NDIM = 4;
typedef DualNumber<RawScalar, NumberVector<NDIM, RawScalar> > FirstDerivType;
typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
typedef SecondDerivType ADType;

using namespace MASA;

template <typename Scalar>
MASA::navierstokes_3d_transient_sutherland<Scalar>::navierstokes_3d_transient_sutherland()
{
  this->mmsname = "navierstokes_3d_transient_sutherland";
  this->dimension = 4;

  this->register_var("L",&L);
  this->register_var("Lt",&Lt);
  this->register_var("rho_0",&rho_0);
  this->register_var("rho_x",&rho_x);
  this->register_var("a_rhox",&a_rhox);
  this->register_var("rho_y",&rho_y);
  this->register_var("a_rhoy",&a_rhoy);
  this->register_var("rho_z",&rho_z);
  this->register_var("a_rhoz",&a_rhoz);
  this->register_var("rho_t",&rho_t);
  this->register_var("a_rhot",&a_rhot);
  this->register_var("p_0",&p_0);
  this->register_var("p_x",&p_x);
  this->register_var("a_px",&a_px);
  this->register_var("p_y",&p_y);
  this->register_var("a_py",&a_py);
  this->register_var("p_z",&p_z);
  this->register_var("a_pz",&a_pz);
  this->register_var("p_t",&p_t);
  this->register_var("a_pt",&a_pt);
  this->register_var("u_0",&u_0);
  this->register_var("u_x",&u_x);
  this->register_var("a_ux",&a_ux);
  this->register_var("u_y",&u_y);
  this->register_var("a_uy",&a_uy);
  this->register_var("u_z",&u_z);
  this->register_var("a_uz",&a_uz);
  this->register_var("u_t",&u_t);
  this->register_var("a_ut",&a_ut);
  this->register_var("v_0",&v_0);
  this->register_var("v_x",&v_x);
  this->register_var("a_vx",&a_vx);
  this->register_var("v_y",&v_y);
  this->register_var("a_vy",&a_vy);
  this->register_var("v_z",&v_z);
  this->register_var("a_vz",&a_vz);
  this->register_var("v_t",&v_t);
  this->register_var("a_vt",&a_vt);
  this->register_var("w_0",&w_0);
  this->register_var("w_x",&w_x);
  this->register_var("a_wx",&a_wx);
  this->register_var("w_y",&w_y);
  this->register_var("a_wy",&a_wy);
  this->register_var("w_z",&w_z);
  this->register_var("a_wz",&a_wz);
  this->register_var("w_t",&w_t);
  this->register_var("a_wt",&a_wt);
  this->register_var("B_mu",&B_mu);
  this->register_var("A_mu",&A_mu);
  this->register_var("Gamma",&Gamma);
  this->register_var("R",&R);
  this->register_var("Pr",&Pr);

  this->init_var();

} // done with constructor

template <typename Scalar>
int MASA::navierstokes_3d_transient_sutherland<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("L",1.5);
  err += this->set_var("Lt",3.);
  err += this->set_var("rho_0",23.);
  err += this->set_var("rho_x",3.);
  err += this->set_var("a_rhox",5.);
  err += this->set_var("rho_y",7.);
  err += this->set_var("a_rhoy",11.);
  err += this->set_var("rho_z",13.);
  err += this->set_var("a_rhoz",17.);
  err += this->set_var("rho_t",1.9);
  err += this->set_var("a_rhot",23);
  err += this->set_var("p_0",29);
  err += this->set_var("p_x",31.);
  err += this->set_var("a_px",37.);
  err += this->set_var("p_y",41.);
  err += this->set_var("a_py",43);
  err += this->set_var("p_z",47.);
  err += this->set_var("a_pz",53);
  err += this->set_var("p_t",6.1);
  err += this->set_var("a_pt",67.);
  err += this->set_var("u_0",2.);
  err += this->set_var("u_x",13.);
  err += this->set_var("a_ux",29.);
  err += this->set_var("u_y",37.);
  err += this->set_var("a_uy",-3.);
  err += this->set_var("u_z",-17.);
  err += this->set_var("a_uz",83.);
  err += this->set_var("u_t",-5.);
  err += this->set_var("a_ut",19.);
  err += this->set_var("v_0",3.);
  err += this->set_var("v_x",11.);
  err += this->set_var("a_vx",23.);
  err += this->set_var("v_y",31.);
  err += this->set_var("a_vy",5.);
  err += this->set_var("v_z",19.);
  err += this->set_var("a_vz",-47.);
  err += this->set_var("v_t",2.);
  err += this->set_var("a_vt",13.);
  err += this->set_var("w_0",5.);
  err += this->set_var("w_x",7.);
  err += this->set_var("a_wx",19.);
  err += this->set_var("w_y",29.);
  err += this->set_var("a_wy",-7.);
  err += this->set_var("w_z",23.);
  err += this->set_var("a_wz",11.);
  err += this->set_var("w_t",11.);
  err += this->set_var("a_wt",17.);
  err += this->set_var("B_mu",110.4);
  err += this->set_var("A_mu",1.458e-6);
  err += this->set_var("Gamma",1.4);
  err += this->set_var("R",287);
  err += this->set_var("Pr",0.7);

  return err;

} // done with init_var

// ----------------------------------------
// Source Terms
// ----------------------------------------

template <typename Scalar>
Scalar MASA::navierstokes_3d_transient_sutherland<Scalar>::eval_q_e(Scalar x, Scalar y, Scalar z, Scalar t)
{
 Scalar RHO;
 Scalar P;
 Scalar U;
 Scalar V;
 Scalar W;
 Scalar T;
 Scalar MU;
 Scalar DMu_Dx;
 Scalar DMu_Dy;
 Scalar DMu_Dz;
 Scalar kappa;
 Scalar Q_e;
 Scalar Q_e_convection;
 Scalar Q_e_work_pressure;
 Scalar Q_e_work_viscous;
 Scalar Q_e_conduction;
 Scalar Q_e_time;
 RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_t * sin(a_rhot * PI * t / Lt);
 U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L) + u_t * cos(a_ut * PI * t / Lt);
 V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_t * sin(a_vt * PI * t / Lt);
 W = w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L) + w_t * cos(a_wt * PI * t / Lt);
 P = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L) + p_t * cos(a_pt * PI * t / Lt);
 T = P / RHO / R;
 MU = A_mu * pow(T, Scalar(0.3e1 / 0.2e1)) / (T + B_mu);
 DMu_Dx = a_rhox * PI * rho_x * MU * MU * cos(a_rhox * PI * x / L) / A_mu / L / RHO / sqrt(T) - 0.3e1 / 0.2e1 * a_rhox * PI * rho_x * MU * cos(a_rhox * PI * x / L) / L / RHO + a_px * PI * p_x * MU * MU * sin(a_px * PI * x / L) / A_mu / L / R / RHO * pow(T, Scalar(-0.3e1 / 0.2e1)) - 0.3e1 / 0.2e1 * a_px * PI * p_x * MU * sin(a_px * PI * x / L) / L / R / RHO / T;
 DMu_Dy = -a_rhoy * PI * rho_y * MU * MU * sin(a_rhoy * PI * y / L) / A_mu / L / RHO / sqrt(T) + 0.3e1 / 0.2e1 * a_rhoy * PI * rho_y * MU * sin(a_rhoy * PI * y / L) / L / RHO - a_py * PI * p_y * MU * MU * cos(a_py * PI * y / L) / A_mu / L / R / RHO * pow(T, Scalar(-0.3e1 / 0.2e1)) + 0.3e1 / 0.2e1 * a_py * PI * p_y * MU * cos(a_py * PI * y / L) / L / R / RHO / T;
 DMu_Dz = a_rhoz * PI * rho_z * MU * MU * cos(a_rhoz * PI * z / L) / A_mu / L / RHO / sqrt(T) - 0.3e1 / 0.2e1 * a_rhoz * PI * rho_z * MU * cos(a_rhoz * PI * z / L) / L / RHO + a_pz * p_z * PI * MU * MU * sin(a_pz * PI * z / L) / A_mu / L / R / RHO * pow(T, Scalar(-0.3e1 / 0.2e1)) - 0.3e1 / 0.2e1 * a_pz * p_z * PI * MU * sin(a_pz * PI * z / L) / L / R / RHO / T;
 kappa = Gamma * R * MU / (Gamma - 0.1e1) / Pr;
 Q_e_convection = a_rhox * PI * rho_x * pow(U, Scalar(0.3e1)) * cos(a_rhox * PI * x / L) / L / 0.2e1 + a_rhox * PI * rho_x * U * V * V * cos(a_rhox * PI * x / L) / L / 0.2e1 + a_rhox * PI * rho_x * U * W * W * cos(a_rhox * PI * x / L) / L / 0.2e1 - a_rhoy * PI * rho_y * U * U * V * sin(a_rhoy * PI * y / L) / L / 0.2e1 - a_rhoy * PI * rho_y * pow(V, Scalar(0.3e1)) * sin(a_rhoy * PI * y / L) / L / 0.2e1 - a_rhoy * PI * rho_y * V * W * W * sin(a_rhoy * PI * y / L) / L / 0.2e1 + a_rhoz * PI * rho_z * U * U * W * cos(a_rhoz * PI * z / L) / L / 0.2e1 + a_rhoz * PI * rho_z * V * V * W * cos(a_rhoz * PI * z / L) / L / 0.2e1 + a_rhoz * PI * rho_z * pow(W, Scalar(0.3e1)) * cos(a_rhoz * PI * z / L) / L / 0.2e1 - a_px * PI * p_x * U * sin(a_px * PI * x / L) / (Gamma - 0.1e1) / L + a_py * PI * p_y * V * cos(a_py * PI * y / L) / (Gamma - 0.1e1) / L - a_pz * p_z * PI * W * sin(a_pz * PI * z / L) / (Gamma - 0.1e1) / L + (0.3e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO * U * U / L / 0.2e1 - (a_uy * u_y * sin(a_uy * PI * y / L) + a_vx * v_x * sin(a_vx * PI * x / L)) * PI * RHO * U * V / L + (-a_uz * u_z * sin(a_uz * PI * z / L) + a_wx * w_x * cos(a_wx * PI * x / L)) * PI * RHO * U * W / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.3e1 * a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO * V * V / L / 0.2e1 + (a_vz * v_z * cos(a_vz * PI * z / L) + a_wy * w_y * cos(a_wy * PI * y / L)) * PI * RHO * V * W / L + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) - 0.3e1 * a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO * W * W / L / 0.2e1 + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * PI * P / (Gamma - 0.1e1) / L;
 Q_e_work_pressure = -a_px * PI * p_x * U * sin(a_px * PI * x / L) / L + a_py * PI * p_y * V * cos(a_py * PI * y / L) / L - a_pz * p_z * PI * W * sin(a_pz * PI * z / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * PI * P / L;
 Q_e_conduction = DMu_Dx * a_rhox * PI * rho_x * kappa * P * cos(a_rhox * PI * x / L) / L / R / MU * pow(RHO, Scalar(-0.2e1)) + DMu_Dx * a_px * PI * p_x * kappa * sin(a_px * PI * x / L) / L / R / MU / RHO - DMu_Dy * a_rhoy * PI * rho_y * kappa * P * sin(a_rhoy * PI * y / L) / L / R / MU * pow(RHO, Scalar(-0.2e1)) - DMu_Dy * a_py * PI * p_y * kappa * cos(a_py * PI * y / L) / L / R / MU / RHO + DMu_Dz * a_rhoz * PI * rho_z * kappa * P * cos(a_rhoz * PI * z / L) / L / R / MU * pow(RHO, Scalar(-0.2e1)) + a_pz * DMu_Dz * p_z * PI * kappa * sin(a_pz * PI * z / L) / L / R / MU / RHO + (a_px * a_px * p_x * cos(a_px * PI * x / L) + a_py * a_py * p_y * sin(a_py * PI * y / L) + a_pz * a_pz * p_z * cos(a_pz * PI * z / L)) * PI * PI * kappa * pow(L, Scalar(-0.2e1)) / R / RHO - (a_rhox * a_rhox * rho_x * sin(a_rhox * PI * x / L) + a_rhoy * a_rhoy * rho_y * cos(a_rhoy * PI * y / L) + a_rhoz * a_rhoz * rho_z * sin(a_rhoz * PI * z / L)) * PI * PI * kappa * P * pow(L, Scalar(-0.2e1)) / R * pow(RHO, Scalar(-0.2e1)) - (0.2e1 * a_rhox * a_px * rho_x * p_x * cos(a_rhox * PI * x / L) * sin(a_px * PI * x / L) + 0.2e1 * a_rhoy * a_py * rho_y * p_y * sin(a_rhoy * PI * y / L) * cos(a_py * PI * y / L) + 0.2e1 * a_pz * a_rhoz * p_z * rho_z * cos(a_rhoz * PI * z / L) * sin(a_pz * PI * z / L)) * PI * PI * kappa * pow(L, Scalar(-0.2e1)) / R * pow(RHO, Scalar(-0.2e1)) - (0.2e1 * a_rhox * a_rhox * rho_x * rho_x * pow(cos(a_rhox * PI * x / L), Scalar(0.2e1)) + 0.2e1 * a_rhoy * a_rhoy * rho_y * rho_y * pow(sin(a_rhoy * PI * y / L), Scalar(0.2e1)) + 0.2e1 * a_rhoz * a_rhoz * rho_z * rho_z * pow(cos(a_rhoz * PI * z / L), Scalar(0.2e1))) * PI * PI * kappa * P * pow(L, Scalar(-0.2e1)) / R * pow(RHO, Scalar(-0.3e1));
 Q_e_work_viscous = (0.4e1 * a_ux * a_ux * u_x * sin(a_ux * PI * x / L) + 0.3e1 * a_uy * a_uy * u_y * cos(a_uy * PI * y / L) + 0.3e1 * a_uz * a_uz * u_z * cos(a_uz * PI * z / L)) * PI * PI * MU * U * pow(L, Scalar(-0.2e1)) / 0.3e1 + (0.3e1 * a_vx * a_vx * v_x * cos(a_vx * PI * x / L) + 0.4e1 * a_vy * a_vy * v_y * sin(a_vy * PI * y / L) + 0.3e1 * a_vz * a_vz * v_z * sin(a_vz * PI * z / L)) * PI * PI * MU * V * pow(L, Scalar(-0.2e1)) / 0.3e1 + (0.3e1 * a_wx * a_wx * w_x * sin(a_wx * PI * x / L) + 0.3e1 * a_wy * a_wy * w_y * sin(a_wy * PI * y / L) + 0.4e1 * a_wz * a_wz * w_z * cos(a_wz * PI * z / L)) * PI * PI * MU * W * pow(L, Scalar(-0.2e1)) / 0.3e1 - 0.2e1 / 0.3e1 * (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) - a_vy * v_y * cos(a_vy * PI * y / L) + a_wz * w_z * sin(a_wz * PI * z / L)) * DMu_Dx * PI * U / L + (a_uy * u_y * sin(a_uy * PI * y / L) + a_vx * v_x * sin(a_vx * PI * x / L)) * DMu_Dx * PI * V / L - (-a_uz * u_z * sin(a_uz * PI * z / L) + a_wx * w_x * cos(a_wx * PI * x / L)) * DMu_Dx * PI * W / L + (a_uy * u_y * sin(a_uy * PI * y / L) + a_vx * v_x * sin(a_vx * PI * x / L)) * DMu_Dy * PI * U / L + 0.2e1 / 0.3e1 * (a_ux * u_x * cos(a_ux * PI * x / L) - 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * DMu_Dy * PI * V / L - (a_vz * v_z * cos(a_vz * PI * z / L) + a_wy * w_y * cos(a_wy * PI * y / L)) * DMu_Dy * PI * W / L - (-a_uz * u_z * sin(a_uz * PI * z / L) + a_wx * w_x * cos(a_wx * PI * x / L)) * DMu_Dz * PI * U / L - (a_vz * v_z * cos(a_vz * PI * z / L) + a_wy * w_y * cos(a_wy * PI * y / L)) * DMu_Dz * PI * V / L + 0.2e1 / 0.3e1 * (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) + 0.2e1 * a_wz * w_z * sin(a_wz * PI * z / L)) * DMu_Dz * PI * W / L - (0.4e1 * a_ux * a_ux * u_x * u_x * pow(cos(a_ux * PI * x / L), Scalar(0.2e1)) - 0.4e1 * a_ux * a_vy * u_x * v_y * cos(a_ux * PI * x / L) * cos(a_vy * PI * y / L) + 0.4e1 * a_ux * a_wz * u_x * w_z * cos(a_ux * PI * x / L) * sin(a_wz * PI * z / L) + 0.3e1 * a_uy * a_uy * u_y * u_y * pow(sin(a_uy * PI * y / L), Scalar(0.2e1)) + 0.6e1 * a_uy * a_vx * u_y * v_x * sin(a_uy * PI * y / L) * sin(a_vx * PI * x / L) + 0.3e1 * a_uz * a_uz * u_z * u_z * pow(sin(a_uz * PI * z / L), Scalar(0.2e1)) - 0.6e1 * a_uz * a_wx * u_z * w_x * sin(a_uz * PI * z / L) * cos(a_wx * PI * x / L) + 0.3e1 * a_vx * a_vx * v_x * v_x * pow(sin(a_vx * PI * x / L), Scalar(0.2e1)) + 0.4e1 * a_vy * a_vy * v_y * v_y * pow(cos(a_vy * PI * y / L), Scalar(0.2e1)) + 0.4e1 * a_vy * a_wz * v_y * w_z * cos(a_vy * PI * y / L) * sin(a_wz * PI * z / L) + 0.3e1 * a_vz * a_vz * v_z * v_z * pow(cos(a_vz * PI * z / L), Scalar(0.2e1)) + 0.6e1 * a_vz * a_wy * v_z * w_y * cos(a_vz * PI * z / L) * cos(a_wy * PI * y / L) + 0.3e1 * a_wx * a_wx * w_x * w_x * pow(cos(a_wx * PI * x / L), Scalar(0.2e1)) + 0.3e1 * a_wy * a_wy * w_y * w_y * pow(cos(a_wy * PI * y / L), Scalar(0.2e1)) + 0.4e1 * a_wz * a_wz * w_z * w_z * pow(sin(a_wz * PI * z / L), Scalar(0.2e1))) * PI * PI * MU * pow(L, Scalar(-0.2e1)) / 0.3e1;
 Q_e_time = -a_ut * u_t * PI * RHO * U * sin(a_ut * PI * t / Lt) / Lt + a_vt * v_t * PI * RHO * V * cos(a_vt * PI * t / Lt) / Lt - a_wt * w_t * PI * RHO * W * sin(a_wt * PI * t / Lt) / Lt + a_rhot * rho_t * PI * U * U * cos(a_rhot * PI * t / Lt) / Lt / 0.2e1 + a_rhot * rho_t * PI * V * V * cos(a_rhot * PI * t / Lt) / Lt / 0.2e1 + a_rhot * rho_t * PI * W * W * cos(a_rhot * PI * t / Lt) / Lt / 0.2e1 - a_pt * p_t * PI * sin(a_pt * PI * t / Lt) / (Gamma - 0.1e1) / Lt;
 Q_e = Q_e_convection + Q_e_work_pressure + Q_e_work_viscous + Q_e_conduction + Q_e_time;
 return(Q_e);
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_transient_sutherland<Scalar>::eval_q_rho(Scalar x, Scalar y, Scalar z, Scalar t)
{
 Scalar RHO;
 Scalar U;
 Scalar V;
 Scalar W;
 Scalar Q_rho;
 Scalar Q_rho_convection;
 Scalar Q_rho_time;
 RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_t * sin(a_rhot * PI * t / Lt);
 U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L) + u_t * cos(a_ut * PI * t / Lt);
 V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_t * sin(a_vt * PI * t / Lt);
 W = w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L) + w_t * cos(a_wt * PI * t / Lt);
 Q_rho_convection = a_rhox * PI * rho_x * U * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * sin(a_rhoy * PI * y / L) / L + a_rhoz * PI * rho_z * W * cos(a_rhoz * PI * z / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO / L;
 Q_rho_time = a_rhot * PI * rho_t * cos(a_rhot * PI * t / Lt) / Lt;
 Q_rho = Q_rho_convection + Q_rho_time;
 return(Q_rho);
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_transient_sutherland<Scalar>::eval_q_u(Scalar x, Scalar y, Scalar z, Scalar t)
{
 Scalar RHO;
 Scalar P;
 Scalar U;
 Scalar V;
 Scalar W;
 Scalar T;
 Scalar MU;
 Scalar DMu_Dx;
 Scalar DMu_Dy;
 Scalar DMu_Dz;
 Scalar Q_u;
 Scalar Q_u_convection;
 Scalar Q_u_pressure;
 Scalar Q_u_viscous;
 Scalar Q_u_time;
 RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_t * sin(a_rhot * PI * t / Lt);
 U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L) + u_t * cos(a_ut * PI * t / Lt);
 V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_t * sin(a_vt * PI * t / Lt);
 W = w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L) + w_t * cos(a_wt * PI * t / Lt);
 P = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L) + p_t * cos(a_pt * PI * t / Lt);
 T = P / RHO / R;
 MU = A_mu * pow(T, Scalar(0.3e1 / 0.2e1)) / (T + B_mu);
 DMu_Dx = a_rhox * PI * rho_x * MU * MU * cos(a_rhox * PI * x / L) / A_mu / L / RHO / sqrt(T) - 0.3e1 / 0.2e1 * a_rhox * PI * rho_x * MU * cos(a_rhox * PI * x / L) / L / RHO + a_px * PI * p_x * MU * MU * sin(a_px * PI * x / L) / A_mu / L / R / RHO * pow(T, Scalar(-0.3e1 / 0.2e1)) - 0.3e1 / 0.2e1 * a_px * PI * p_x * MU * sin(a_px * PI * x / L) / L / R / RHO / T;
 DMu_Dy = -a_rhoy * PI * rho_y * MU * MU * sin(a_rhoy * PI * y / L) / A_mu / L / RHO / sqrt(T) + 0.3e1 / 0.2e1 * a_rhoy * PI * rho_y * MU * sin(a_rhoy * PI * y / L) / L / RHO - a_py * PI * p_y * MU * MU * cos(a_py * PI * y / L) / A_mu / L / R / RHO * pow(T, Scalar(-0.3e1 / 0.2e1)) + 0.3e1 / 0.2e1 * a_py * PI * p_y * MU * cos(a_py * PI * y / L) / L / R / RHO / T;
 DMu_Dz = a_rhoz * PI * rho_z * MU * MU * cos(a_rhoz * PI * z / L) / A_mu / L / RHO / sqrt(T) - 0.3e1 / 0.2e1 * a_rhoz * PI * rho_z * MU * cos(a_rhoz * PI * z / L) / L / RHO + a_pz * p_z * PI * MU * MU * sin(a_pz * PI * z / L) / A_mu / L / R / RHO * pow(T, Scalar(-0.3e1 / 0.2e1)) - 0.3e1 / 0.2e1 * a_pz * p_z * PI * MU * sin(a_pz * PI * z / L) / L / R / RHO / T;
 Q_u_convection = a_rhox * PI * rho_x * U * U * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * U * V * sin(a_rhoy * PI * y / L) / L + a_rhoz * PI * rho_z * U * W * cos(a_rhoz * PI * z / L) / L - a_uy * PI * u_y * RHO * V * sin(a_uy * PI * y / L) / L - a_uz * PI * u_z * RHO * W * sin(a_uz * PI * z / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO * U / L;
 Q_u_pressure = -a_px * PI * p_x * sin(a_px * PI * x / L) / L;
 Q_u_viscous = (0.4e1 * a_ux * a_ux * u_x * sin(a_ux * PI * x / L) + 0.3e1 * a_uy * a_uy * u_y * cos(a_uy * PI * y / L) + 0.3e1 * a_uz * a_uz * u_z * cos(a_uz * PI * z / L)) * PI * PI * MU * pow(L, Scalar(-0.2e1)) / 0.3e1 - 0.2e1 / 0.3e1 * (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) - a_vy * v_y * cos(a_vy * PI * y / L) + a_wz * w_z * sin(a_wz * PI * z / L)) * DMu_Dx * PI / L + (a_uy * u_y * sin(a_uy * PI * y / L) + a_vx * v_x * sin(a_vx * PI * x / L)) * DMu_Dy * PI / L + (a_uz * u_z * sin(a_uz * PI * z / L) - a_wx * w_x * cos(a_wx * PI * x / L)) * DMu_Dz * PI / L;
 Q_u_time = -a_ut * u_t * PI * RHO * sin(a_ut * PI * t / Lt) / Lt + a_rhot * rho_t * PI * U * cos(a_rhot * PI * t / Lt) / Lt;
 Q_u = Q_u_convection + Q_u_pressure + Q_u_viscous + Q_u_time;
 return(Q_u);
}
template <typename Scalar>
Scalar MASA::navierstokes_3d_transient_sutherland<Scalar>::eval_q_v(Scalar x, Scalar y, Scalar z, Scalar t)
{
 Scalar RHO;
 Scalar P;
 Scalar U;
 Scalar V;
 Scalar W;
 Scalar T;
 Scalar MU;
 Scalar DMu_Dx;
 Scalar DMu_Dy;
 Scalar DMu_Dz;
 Scalar Q_v;
 Scalar Q_v_convection;
 Scalar Q_v_pressure;
 Scalar Q_v_viscous;
 Scalar Q_v_time;
 RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_t * sin(a_rhot * PI * t / Lt);
 U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L) + u_t * cos(a_ut * PI * t / Lt);
 V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_t * sin(a_vt * PI * t / Lt);
 W = w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L) + w_t * cos(a_wt * PI * t / Lt);
 P = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L) + p_t * cos(a_pt * PI * t / Lt);
 T = P / RHO / R;
 MU = A_mu * pow(T, Scalar(0.3e1 / 0.2e1)) / (T + B_mu);
 DMu_Dx = a_rhox * PI * rho_x * MU * MU * cos(a_rhox * PI * x / L) / A_mu / L / RHO / sqrt(T) - 0.3e1 / 0.2e1 * a_rhox * PI * rho_x * MU * cos(a_rhox * PI * x / L) / L / RHO + a_px * PI * p_x * MU * MU * sin(a_px * PI * x / L) / A_mu / L / R / RHO * pow(T, Scalar(-0.3e1 / 0.2e1)) - 0.3e1 / 0.2e1 * a_px * PI * p_x * MU * sin(a_px * PI * x / L) / L / R / RHO / T;
 DMu_Dy = -a_rhoy * PI * rho_y * MU * MU * sin(a_rhoy * PI * y / L) / A_mu / L / RHO / sqrt(T) + 0.3e1 / 0.2e1 * a_rhoy * PI * rho_y * MU * sin(a_rhoy * PI * y / L) / L / RHO - a_py * PI * p_y * MU * MU * cos(a_py * PI * y / L) / A_mu / L / R / RHO * pow(T, Scalar(-0.3e1 / 0.2e1)) + 0.3e1 / 0.2e1 * a_py * PI * p_y * MU * cos(a_py * PI * y / L) / L / R / RHO / T;
 DMu_Dz = a_rhoz * PI * rho_z * MU * MU * cos(a_rhoz * PI * z / L) / A_mu / L / RHO / sqrt(T) - 0.3e1 / 0.2e1 * a_rhoz * PI * rho_z * MU * cos(a_rhoz * PI * z / L) / L / RHO + a_pz * p_z * PI * MU * MU * sin(a_pz * PI * z / L) / A_mu / L / R / RHO * pow(T, Scalar(-0.3e1 / 0.2e1)) - 0.3e1 / 0.2e1 * a_pz * p_z * PI * MU * sin(a_pz * PI * z / L) / L / R / RHO / T;
 Q_v_convection = a_rhox * PI * rho_x * U * V * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * V * sin(a_rhoy * PI * y / L) / L + a_rhoz * PI * rho_z * V * W * cos(a_rhoz * PI * z / L) / L - a_vx * PI * v_x * RHO * U * sin(a_vx * PI * x / L) / L + a_vz * PI * v_z * RHO * W * cos(a_vz * PI * z / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO * V / L;
 Q_v_pressure = a_py * PI * p_y * cos(a_py * PI * y / L) / L;
 Q_v_viscous = (0.3e1 * a_vx * a_vx * v_x * cos(a_vx * PI * x / L) + 0.4e1 * a_vy * a_vy * v_y * sin(a_vy * PI * y / L) + 0.3e1 * a_vz * a_vz * v_z * sin(a_vz * PI * z / L)) * PI * PI * MU * pow(L, Scalar(-0.2e1)) / 0.3e1 + (a_uy * u_y * sin(a_uy * PI * y / L) + a_vx * v_x * sin(a_vx * PI * x / L)) * DMu_Dx * PI / L + 0.2e1 / 0.3e1 * (a_ux * u_x * cos(a_ux * PI * x / L) - 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * DMu_Dy * PI / L - (a_vz * v_z * cos(a_vz * PI * z / L) + a_wy * w_y * cos(a_wy * PI * y / L)) * DMu_Dz * PI / L;
 Q_v_time = a_vt * v_t * PI * RHO * cos(a_vt * PI * t / Lt) / Lt + a_rhot * rho_t * PI * V * cos(a_rhot * PI * t / Lt) / Lt;
 Q_v = Q_v_convection + Q_v_pressure + Q_v_viscous + Q_v_time;
 return(Q_v);
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_transient_sutherland<Scalar>::eval_q_w(Scalar x, Scalar y, Scalar z, Scalar t)
{
 Scalar RHO;
 Scalar P;
 Scalar U;
 Scalar V;
 Scalar W;
 Scalar T;
 Scalar MU;
 Scalar DMu_Dx;
 Scalar DMu_Dy;
 Scalar DMu_Dz;
 Scalar Q_w;
 Scalar Q_w_convection;
 Scalar Q_w_pressure;
 Scalar Q_w_viscous;
 Scalar Q_w_time;
 RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_t * sin(a_rhot * PI * t / Lt);
 U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L) + u_t * cos(a_ut * PI * t / Lt);
 V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_t * sin(a_vt * PI * t / Lt);
 W = w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L) + w_t * cos(a_wt * PI * t / Lt);
 P = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L) + p_t * cos(a_pt * PI * t / Lt);
 T = P / RHO / R;
 MU = A_mu * pow(T, Scalar(0.3e1 / 0.2e1)) / (T + B_mu);
 DMu_Dx = a_rhox * PI * rho_x * MU * MU * cos(a_rhox * PI * x / L) / A_mu / L / RHO / sqrt(T) - 0.3e1 / 0.2e1 * a_rhox * PI * rho_x * MU * cos(a_rhox * PI * x / L) / L / RHO + a_px * PI * p_x * MU * MU * sin(a_px * PI * x / L) / A_mu / L / R / RHO * pow(T, Scalar(-0.3e1 / 0.2e1)) - 0.3e1 / 0.2e1 * a_px * PI * p_x * MU * sin(a_px * PI * x / L) / L / R / RHO / T;
 DMu_Dy = -a_rhoy * PI * rho_y * MU * MU * sin(a_rhoy * PI * y / L) / A_mu / L / RHO / sqrt(T) + 0.3e1 / 0.2e1 * a_rhoy * PI * rho_y * MU * sin(a_rhoy * PI * y / L) / L / RHO - a_py * PI * p_y * MU * MU * cos(a_py * PI * y / L) / A_mu / L / R / RHO * pow(T, Scalar(-0.3e1 / 0.2e1)) + 0.3e1 / 0.2e1 * a_py * PI * p_y * MU * cos(a_py * PI * y / L) / L / R / RHO / T;
 DMu_Dz = a_rhoz * PI * rho_z * MU * MU * cos(a_rhoz * PI * z / L) / A_mu / L / RHO / sqrt(T) - 0.3e1 / 0.2e1 * a_rhoz * PI * rho_z * MU * cos(a_rhoz * PI * z / L) / L / RHO + a_pz * p_z * PI * MU * MU * sin(a_pz * PI * z / L) / A_mu / L / R / RHO * pow(T, Scalar(-0.3e1 / 0.2e1)) - 0.3e1 / 0.2e1 * a_pz * p_z * PI * MU * sin(a_pz * PI * z / L) / L / R / RHO / T;
 Q_w_convection = a_rhox * PI * rho_x * U * W * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * W * sin(a_rhoy * PI * y / L) / L + a_rhoz * PI * rho_z * W * W * cos(a_rhoz * PI * z / L) / L + a_wx * PI * w_x * RHO * U * cos(a_wx * PI * x / L) / L + a_wy * PI * w_y * RHO * V * cos(a_wy * PI * y / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) - 0.2e1 * a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO * W / L;
 Q_w_pressure = -a_pz * p_z * PI * sin(a_pz * PI * z / L) / L;
 Q_w_viscous = (0.3e1 * a_wx * a_wx * w_x * sin(a_wx * PI * x / L) + 0.3e1 * a_wy * a_wy * w_y * sin(a_wy * PI * y / L) + 0.4e1 * a_wz * a_wz * w_z * cos(a_wz * PI * z / L)) * PI * PI * MU * pow(L, Scalar(-0.2e1)) / 0.3e1 - (-a_uz * u_z * sin(a_uz * PI * z / L) + a_wx * w_x * cos(a_wx * PI * x / L)) * DMu_Dx * PI / L - (a_vz * v_z * cos(a_vz * PI * z / L) + a_wy * w_y * cos(a_wy * PI * y / L)) * DMu_Dy * PI / L + 0.2e1 / 0.3e1 * (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) + 0.2e1 * a_wz * w_z * sin(a_wz * PI * z / L)) * DMu_Dz * PI / L;
 Q_w_time = -a_wt * w_t * PI * RHO * sin(a_wt * PI * t / Lt) / Lt + a_rhot * rho_t * PI * W * cos(a_rhot * PI * t / Lt) / Lt;
 Q_w = Q_w_convection + Q_w_pressure + Q_w_viscous + Q_w_time;
 return(Q_w);
}


// ----------------------------------------
// Analytical Terms
// ----------------------------------------
template <typename Scalar>
Scalar MASA::navierstokes_3d_transient_sutherland<Scalar>::eval_exact_rho(Scalar x, Scalar y, Scalar z, Scalar t)
{
  using std::cos;
  using std::sin;
  
  Scalar rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot * pi * t / Lt);
  
  return rho_an;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_transient_sutherland<Scalar>::eval_exact_p(Scalar x, Scalar y, Scalar z, Scalar t)
{
  using std::cos;
  using std::sin;
  
  Scalar p_an = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt);
  
  return p_an;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_transient_sutherland<Scalar>::eval_exact_u(Scalar x, Scalar y, Scalar z, Scalar t)
{
  using std::cos;
  using std::sin;
  
  Scalar u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt);
  
  return u_an;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_transient_sutherland<Scalar>::eval_exact_v(Scalar x, Scalar y, Scalar z, Scalar t)
{
  using std::cos;
  using std::sin;
  
  Scalar v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt);
  
  return v_an;
}


template <typename Scalar>
Scalar MASA::navierstokes_3d_transient_sutherland<Scalar>::eval_exact_w(Scalar x, Scalar y, Scalar z, Scalar t)
{
  using std::cos;
  using std::sin;
  
  Scalar w_an = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt);
  
  return w_an;
}






// ----------------------------------------
// Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::navierstokes_3d_transient_sutherland);


#endif // HAVE_METAPHYSICL
