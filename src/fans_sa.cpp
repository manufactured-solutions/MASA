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
// $Id: heat.cpp 18162 2011-03-01 05:23:07Z nick $
//
// fans_sa.cpp: these are class definitions and methods of favre averaged
//              navier stokes with spelart alamaras 
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h> 

using namespace MASA;

/* ------------------------------------------------
 *
 *   favre averaged navier stokes transient
 *            infinite distance to wall
 *
 * -----------------------------------------------
 */ 
template <typename Scalar>
MASA::fans_sa_transient_free_shear<Scalar>::fans_sa_transient_free_shear()
{
  this->mmsname = "fans_sa_transient_free_shear";
  this->dimension=2;

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
  this->register_var("mu",&mu);

  this->register_var("u_t",&u_t);
  this->register_var("v_t",&v_t);
  this->register_var("p_t",&p_t);
  this->register_var("rho_t",&rho_t);
  this->register_var("a_ut",&a_ut);
  this->register_var("a_vt",&a_vt);
  this->register_var("a_pt",&a_pt);
  this->register_var("a_rhot",&a_rhot);

  this->register_var("nu_sa_0",&nu_sa_0);
  this->register_var("nu_sa_x",&nu_sa_x);
  this->register_var("nu_sa_y",&nu_sa_y);
  this->register_var("nu_sa_t",&nu_sa_t);

  this->register_var("a_nusax",&a_nusax);
  this->register_var("a_nusay",&a_nusay);
  this->register_var("a_nusat",&a_nusat);

  this->register_var("c_v1",&c_v1);
  this->register_var("c_b1",&c_b1);
  this->register_var("c_b2",&c_b2);
  this->register_var("c_w1",&c_w1);
  this->register_var("c_w2",&c_w2);
  this->register_var("c_w3",&c_w3);
  this->register_var("kappa",&kappa);
  this->register_var("sigma",&sigma);
  //this->register_var("cp",&cp);
  //this->register_var("cv",&cv);
  this->register_var("Gamma", &Gamma);
  this->register_var("R", &R);
  this->register_var("Pr",&Pr);
  this->register_var("Pr_t",&Pr_t);

  this->init_var();

}

template <typename Scalar>
int MASA::fans_sa_transient_free_shear<Scalar>::init_var()
{
  int err = 0;
  
  err += this->set_var("u_0",10.0);
  err += this->set_var("u_x",1.0);
  err += this->set_var("u_y",8.0);
  err += this->set_var("v_0",0.0);
  err += this->set_var("v_x",0.0);
  err += this->set_var("v_y",1.0);
  err += this->set_var("rho_0",1.0);
  err += this->set_var("rho_x",0.1);
  err += this->set_var("rho_y",-0.2);
  err += this->set_var("p_0",1.0e5); //this->set_var("p_0",1.0);
  err += this->set_var("p_x",10.0); //this->set_var("p_x",0.1);
  err += this->set_var("p_y",10.0); //this->set_var("p_y",0.1);
  err += this->set_var("a_px",2.0);
  err += this->set_var("a_py",1.0);
  err += this->set_var("a_rhox",1.0);
  err += this->set_var("a_rhoy",1.0);
  err += this->set_var("a_ux",3.0);
  err += this->set_var("a_uy",1.0);
  err += this->set_var("a_vx",2.0);
  err += this->set_var("a_vy",0.5);
  err += this->set_var("mu",1.0e-3);
  err += this->set_var("L",1.0);

  err += this->set_var("nu_sa_0",0.2);
  err += this->set_var("nu_sa_x",0.1);
  err += this->set_var("nu_sa_y",0.2);
  err += this->set_var("nu_sa_t",0.0);

  err += this->set_var("a_nusay",0.5);
  err += this->set_var("a_nusax",1.0);
  err += this->set_var("a_nusat",0.0);

  err += this->set_var("c_v1",7.1);
  err += this->set_var("c_b1",0.1355);
  err += this->set_var("c_b2",0.622);
  err += this->set_var("c_w1",3.23906781677573e+00); // = c_b1/(kappa*kappa) + (1.0+c_b2)/sigma;
  err += this->set_var("c_w2",0.3);
  err += this->set_var("c_w3",2.0);
  err += this->set_var("kappa",0.41);
  err += this->set_var("sigma",2./3.);
  //err += this->set_var("cp",3.5*8314.472/28.96); // = 3.5*Runiversal/(mol wght for air)
  //err += this->set_var("cv",3.5*8314.472/28.96/1.4); // = cp/gamma = cp/1.4
  err += this->set_var("Gamma", 1.4);
  err += this->set_var("R", 287.0);
  
  err += this->set_var("Pr",0.71);
  err += this->set_var("Pr_t",0.9);

  err += this->set_var("u_t",0.0);
  err += this->set_var("v_t",0.0);
  err += this->set_var("p_t",0.0);
  err += this->set_var("rho_t",0.0);
  err += this->set_var("a_ut",0.0);
  err += this->set_var("a_vt",0.0);
  err += this->set_var("a_pt",0.0);
  err += this->set_var("a_rhot",0.0);
  
  return err;  
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_free_shear<Scalar>::eval_q_rho_u(Scalar x,Scalar y)
{
  return eval_q_rho_u(x,y,0.0);
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_free_shear<Scalar>::eval_q_rho_u(Scalar x,Scalar y,Scalar t)
{
  Scalar Q_u;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar NU_SA;
  Scalar f_v1;
  Scalar chi;
  NU_SA = nu_sa_0 + nu_sa_x * cos(a_nusax * PI * x / L) + nu_sa_y * cos(a_nusay * PI * y / L) + nu_sa_t * cos(a_nusat * PI * t / L);
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  chi = RHO * NU_SA / mu;
  f_v1 = pow(chi, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1));
  Q_u = a_rhox * PI * rho_x * U * U * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * U * V * sin(a_rhoy * PI * y / L) / L - a_uy * PI * u_y * RHO * V * sin(a_uy * PI * y / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * U / L - a_px * PI * p_x * sin(a_px * PI * x / L) / L + (0.4e1 / 0.3e1 * a_ux * a_ux * u_x * sin(a_ux * PI * x / L) + a_uy * a_uy * u_y * cos(a_uy * PI * y / L)) * f_v1 * PI * PI * RHO * NU_SA * pow(L, -0.2e1) + (0.4e1 / 0.3e1 * a_ux * a_nusax * u_x * nu_sa_x * cos(a_ux * PI * x / L) * sin(a_nusax * PI * x / L) - a_uy * a_nusay * u_y * nu_sa_y * sin(a_uy * PI * y / L) * sin(a_nusay * PI * y / L) - a_vx * a_nusay * v_x * nu_sa_y * sin(a_vx * PI * x / L) * sin(a_nusay * PI * y / L) - 0.2e1 / 0.3e1 * a_vy * a_nusax * v_y * nu_sa_x * cos(a_vy * PI * y / L) * sin(a_nusax * PI * x / L)) * f_v1 * PI * PI * RHO * pow(L, -0.2e1) + (-0.4e1 / 0.3e1 * a_rhox * a_ux * rho_x * u_x * cos(a_rhox * PI * x / L) * cos(a_ux * PI * x / L) + 0.2e1 / 0.3e1 * a_rhox * a_vy * rho_x * v_y * cos(a_rhox * PI * x / L) * cos(a_vy * PI * y / L) - a_rhoy * a_uy * rho_y * u_y * sin(a_rhoy * PI * y / L) * sin(a_uy * PI * y / L) - a_rhoy * a_vx * rho_y * v_x * sin(a_rhoy * PI * y / L) * sin(a_vx * PI * x / L)) * f_v1 * PI * PI * NU_SA * pow(L, -0.2e1) + (pow(c_v1, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1)) + f_v1) * (0.4e1 * a_ux * a_ux * u_x * sin(a_ux * PI * x / L) + 0.3e1 * a_uy * a_uy * u_y * cos(a_uy * PI * y / L)) * PI * PI * mu * pow(L, -0.2e1) / 0.3e1 + a_rhot * PI * rho_t * U * cos(a_rhot * PI * t / L) / L - a_ut * PI * u_t * RHO * sin(a_ut * PI * t / L) / L;
  return(Q_u);

}

template <typename Scalar>
Scalar MASA::fans_sa_transient_free_shear<Scalar>::eval_q_rho_v(Scalar x,Scalar y)
{
  return eval_q_rho_v(x,y,0.0);
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_free_shear<Scalar>::eval_q_rho_v(Scalar x,Scalar y,Scalar t)
{
  Scalar Q_v;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar NU_SA;
  Scalar f_v1;
  Scalar chi;
  NU_SA = nu_sa_0 + nu_sa_x * cos(a_nusax * PI * x / L) + nu_sa_y * cos(a_nusay * PI * y / L) + nu_sa_t * cos(a_nusat * PI * t / L);
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  chi = RHO * NU_SA / mu;
  f_v1 = pow(chi, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1));
  Q_v = a_rhox * PI * rho_x * U * V * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * V * sin(a_rhoy * PI * y / L) / L - a_vx * PI * v_x * RHO * U * sin(a_vx * PI * x / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * V / L + a_py * PI * p_y * cos(a_py * PI * y / L) / L + (a_vx * a_vx * v_x * cos(a_vx * PI * x / L) + 0.4e1 / 0.3e1 * a_vy * a_vy * v_y * sin(a_vy * PI * y / L)) * f_v1 * PI * PI * RHO * NU_SA * pow(L, -0.2e1) + (-0.2e1 / 0.3e1 * a_ux * a_nusay * u_x * nu_sa_y * cos(a_ux * PI * x / L) * sin(a_nusay * PI * y / L) - a_uy * a_nusax * u_y * nu_sa_x * sin(a_uy * PI * y / L) * sin(a_nusax * PI * x / L) - a_vx * a_nusax * v_x * nu_sa_x * sin(a_vx * PI * x / L) * sin(a_nusax * PI * x / L) + 0.4e1 / 0.3e1 * a_vy * a_nusay * v_y * nu_sa_y * cos(a_vy * PI * y / L) * sin(a_nusay * PI * y / L)) * f_v1 * PI * PI * RHO * pow(L, -0.2e1) + (a_rhox * a_uy * rho_x * u_y * cos(a_rhox * PI * x / L) * sin(a_uy * PI * y / L) + a_rhox * a_vx * rho_x * v_x * cos(a_rhox * PI * x / L) * sin(a_vx * PI * x / L) - 0.2e1 / 0.3e1 * a_rhoy * a_ux * rho_y * u_x * sin(a_rhoy * PI * y / L) * cos(a_ux * PI * x / L) + 0.4e1 / 0.3e1 * a_rhoy * a_vy * rho_y * v_y * sin(a_rhoy * PI * y / L) * cos(a_vy * PI * y / L)) * f_v1 * PI * PI * NU_SA * pow(L, -0.2e1) + (pow(c_v1, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1)) + f_v1) * (0.3e1 * a_vx * a_vx * v_x * cos(a_vx * PI * x / L) + 0.4e1 * a_vy * a_vy * v_y * sin(a_vy * PI * y / L)) * PI * PI * mu * pow(L, -0.2e1) / 0.3e1 + a_rhot * PI * rho_t * V * cos(a_rhot * PI * t / L) / L + a_vt * PI * v_t * RHO * cos(a_vt * PI * t / L) / L;
  return(Q_v);
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_free_shear<Scalar>::eval_q_rho(Scalar x,Scalar y)
{
  return eval_q_rho(x,y,0.0);
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_free_shear<Scalar>::eval_q_rho(Scalar x,Scalar y,Scalar t)
{
  Scalar Q_rho;
  Scalar RHO;
  Scalar U;
  Scalar V;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Q_rho = a_rhox * PI * rho_x * U * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * sin(a_rhoy * PI * y / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO / L + a_rhot * PI * rho_t * cos(a_rhot * PI * t / L) / L;
  return(Q_rho);
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_free_shear<Scalar>::eval_q_nu(Scalar x,Scalar y)
{
  return eval_q_nu(x,y,0.0);
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_free_shear<Scalar>::eval_q_nu(Scalar x,Scalar y,Scalar t)
{
  Scalar Q_nu;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar NU_SA;
  NU_SA = nu_sa_0 + nu_sa_x * cos(a_nusax * PI * x / L) + nu_sa_y * cos(a_nusay * PI * y / L) + nu_sa_t * cos(a_nusat * PI * t / L);
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Q_nu = a_rhox * PI * rho_x * U * NU_SA * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * NU_SA * sin(a_rhoy * PI * y / L) / L - a_nusax * PI * nu_sa_x * RHO * U * sin(a_nusax * PI * x / L) / L - a_nusay * PI * nu_sa_y * RHO * V * sin(a_nusay * PI * y / L) / L + a_rhot * PI * rho_t * NU_SA * cos(a_rhot * PI * t / L) / L - a_nusat * PI * nu_sa_t * RHO * sin(a_nusat * PI * t / L) / L - (a_nusax * a_nusax * nu_sa_x * nu_sa_x * pow(sin(a_nusax * PI * x / L), 0.2e1) + a_nusay * a_nusay * nu_sa_y * nu_sa_y * pow(sin(a_nusay * PI * y / L), 0.2e1)) * (0.1e1 + c_b2) * PI * PI * RHO / sigma * pow(L, -0.2e1) + (a_rhox * a_nusax * rho_x * nu_sa_x * cos(a_rhox * PI * x / L) * sin(a_nusax * PI * x / L) - a_rhoy * a_nusay * rho_y * nu_sa_y * sin(a_rhoy * PI * y / L) * sin(a_nusay * PI * y / L)) * PI * PI * NU_SA / sigma * pow(L, -0.2e1) - c_b1 * sqrt(pow(-a_uy * u_y * sin(a_uy * PI * y / L) + a_vx * v_x * sin(a_vx * PI * x / L), 0.2e1) * pow(L, -0.2e1)) * PI * RHO * NU_SA + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * NU_SA / L + (a_nusax * a_nusax * nu_sa_x * cos(a_nusax * PI * x / L) + a_nusay * a_nusay * nu_sa_y * cos(a_nusay * PI * y / L)) * (RHO * NU_SA + mu) * PI * PI / sigma * pow(L, -0.2e1);
  return(Q_nu);
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_free_shear<Scalar>::eval_q_rho_e(Scalar x,Scalar y)
{
  return eval_q_rho_e(x,y,0.0);
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_free_shear<Scalar>::eval_q_rho_e(Scalar x,Scalar y,Scalar t)
{
  double Q_E;
  double RHO;
  double U;
  double V;
  double P;
  double NU_SA;
  double chi;
  double f_v1;
  //double R;
  double cp;
  double cv;
  double mu_t;
  NU_SA = nu_sa_0 + nu_sa_x * cos(a_nusax * PI * x / L) + nu_sa_y * cos(a_nusay * PI * y / L) + nu_sa_t * cos(a_nusat * PI * t / L);
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  P = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_t * cos(a_pt * PI * t / L);
  chi = RHO * NU_SA / mu;
  f_v1 = pow(chi, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1));
  mu_t = RHO * NU_SA * f_v1;
  //R = cp - cv;
  cv = R/(Gamma-1.0);
  cp = Gamma*cv;
  Q_E = -(mu_t / Pr_t + mu / Pr) * (-(a_px * a_px * p_x * cos(a_px * PI * x / L) + a_py * a_py * p_y * sin(a_py * PI * y / L)) * cp * PI * PI * pow(L, -0.2e1) / R / RHO + (a_rhox * a_rhox * rho_x * sin(a_rhox * PI * x / L) + a_rhoy * a_rhoy * rho_y * cos(a_rhoy * PI * y / L)) * cp * PI * PI * P * pow(L, -0.2e1) / R * pow(RHO, -0.2e1)) + (U * U + V * V) * a_rhox * PI * rho_x * U * cos(a_rhox * PI * x / L) / L / 0.2e1 - (a_rhox * a_rhox * rho_x * rho_x * pow(cos(a_rhox * PI * x / L), 0.2e1) + a_rhoy * a_rhoy * rho_y * rho_y * pow(sin(a_rhoy * PI * y / L), 0.2e1)) * (Pr * mu_t + 0.2e1 * Pr_t * mu) * cp * PI * PI * P / Pr / Pr_t * pow(L, -0.2e1) / R * pow(RHO, -0.3e1) - (a_rhox * a_px * rho_x * p_x * cos(a_rhox * PI * x / L) * sin(a_px * PI * x / L) + a_rhoy * a_py * rho_y * p_y * sin(a_rhoy * PI * y / L) * cos(a_py * PI * y / L)) * (Pr * mu_t + 0.2e1 * Pr_t * mu) * cp * PI * PI / Pr / Pr_t * pow(L, -0.2e1) / R * pow(RHO, -0.2e1) - (a_rhox * a_nusax * rho_x * nu_sa_x * cos(a_rhox * PI * x / L) * sin(a_nusax * PI * x / L) - a_rhoy * a_nusay * rho_y * nu_sa_y * sin(a_rhoy * PI * y / L) * sin(a_nusay * PI * y / L)) * cp * PI * PI * mu_t * P * pow(L, -0.2e1) / Pr_t / R * pow(RHO, -0.2e1) / NU_SA - (a_px * a_nusax * p_x * nu_sa_x * sin(a_px * PI * x / L) * sin(a_nusax * PI * x / L) - a_py * a_nusay * p_y * nu_sa_y * cos(a_py * PI * y / L) * sin(a_nusay * PI * y / L)) * cp * PI * PI * mu_t * pow(L, -0.2e1) / Pr_t / R / RHO / NU_SA + (0.4e1 / 0.3e1 * a_ux * a_nusax * u_x * nu_sa_x * cos(a_ux * PI * x / L) * sin(a_nusax * PI * x / L) - a_uy * a_nusay * u_y * nu_sa_y * sin(a_uy * PI * y / L) * sin(a_nusay * PI * y / L) - a_vx * a_nusay * v_x * nu_sa_y * sin(a_vx * PI * x / L) * sin(a_nusay * PI * y / L) - 0.2e1 / 0.3e1 * a_vy * a_nusax * v_y * nu_sa_x * cos(a_vy * PI * y / L) * sin(a_nusax * PI * x / L)) * PI * PI * f_v1 * RHO * U * pow(L, -0.2e1) + (-0.2e1 / 0.3e1 * a_ux * a_nusay * u_x * nu_sa_y * cos(a_ux * PI * x / L) * sin(a_nusay * PI * y / L) - a_uy * a_nusax * u_y * nu_sa_x * sin(a_uy * PI * y / L) * sin(a_nusax * PI * x / L) - a_vx * a_nusax * v_x * nu_sa_x * sin(a_vx * PI * x / L) * sin(a_nusax * PI * x / L) + 0.4e1 / 0.3e1 * a_vy * a_nusay * v_y * nu_sa_y * cos(a_vy * PI * y / L) * sin(a_nusay * PI * y / L)) * PI * PI * f_v1 * RHO * V * pow(L, -0.2e1) + (-0.4e1 / 0.3e1 * a_rhox * a_ux * rho_x * u_x * cos(a_rhox * PI * x / L) * cos(a_ux * PI * x / L) + 0.2e1 / 0.3e1 * a_rhox * a_vy * rho_x * v_y * cos(a_rhox * PI * x / L) * cos(a_vy * PI * y / L) - a_rhoy * a_uy * rho_y * u_y * sin(a_rhoy * PI * y / L) * sin(a_uy * PI * y / L) - a_rhoy * a_vx * rho_y * v_x * sin(a_rhoy * PI * y / L) * sin(a_vx * PI * x / L)) * PI * PI * f_v1 * U * NU_SA * pow(L, -0.2e1) + (a_rhox * a_uy * rho_x * u_y * cos(a_rhox * PI * x / L) * sin(a_uy * PI * y / L) + a_rhox * a_vx * rho_x * v_x * cos(a_rhox * PI * x / L) * sin(a_vx * PI * x / L) - 0.2e1 / 0.3e1 * a_rhoy * a_ux * rho_y * u_x * sin(a_rhoy * PI * y / L) * cos(a_ux * PI * x / L) + 0.4e1 / 0.3e1 * a_rhoy * a_vy * rho_y * v_y * sin(a_rhoy * PI * y / L) * cos(a_vy * PI * y / L)) * PI * PI * f_v1 * V * NU_SA * pow(L, -0.2e1) + (U * U + V * V) * a_rhot * PI * rho_t * cos(a_rhot * PI * t / L) / L / 0.2e1 + (f_v1 + pow(c_v1, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1))) * (0.4e1 * a_ux * a_ux * u_x * sin(a_ux * PI * x / L) + 0.3e1 * a_uy * a_uy * u_y * cos(a_uy * PI * y / L)) * PI * PI * mu * U * pow(L, -0.2e1) / 0.3e1 + (f_v1 + pow(c_v1, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1))) * (0.3e1 * a_vx * a_vx * v_x * cos(a_vx * PI * x / L) + 0.4e1 * a_vy * a_vy * v_y * sin(a_vy * PI * y / L)) * PI * PI * mu * V * pow(L, -0.2e1) / 0.3e1 - (a_uy * u_y * sin(a_uy * PI * y / L) + a_vx * v_x * sin(a_vx * PI * x / L)) * PI * RHO * U * V / L + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * cp * PI * P / L / R - a_ut * PI * u_t * RHO * U * sin(a_ut * PI * t / L) / L - (U * U + V * V) * a_rhoy * PI * rho_y * V * sin(a_rhoy * PI * y / L) / L / 0.2e1 + a_vt * PI * v_t * RHO * V * cos(a_vt * PI * t / L) / L + cp * a_py * PI * p_y * V * cos(a_py * PI * y / L) / L / R - cp * a_px * PI * p_x * U * sin(a_px * PI * x / L) / L / R + (0.4e1 * a_ux * a_ux * u_x * sin(a_ux * PI * x / L) + 0.3e1 * a_uy * a_uy * u_y * cos(a_uy * PI * y / L)) * PI * PI * mu_t * U * pow(L, -0.2e1) / 0.3e1 + (0.3e1 * a_vx * a_vx * v_x * cos(a_vx * PI * x / L) + 0.4e1 * a_vy * a_vy * v_y * sin(a_vy * PI * y / L)) * PI * PI * mu_t * V * pow(L, -0.2e1) / 0.3e1 + (0.3e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * U * U / L / 0.2e1 + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.3e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * V * V / L / 0.2e1 - (f_v1 + pow(c_v1, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1))) * (0.4e1 * a_ux * a_ux * u_x * u_x * pow(cos(a_ux * PI * x / L), 0.2e1) - 0.4e1 * a_ux * a_vy * u_x * v_y * cos(a_ux * PI * x / L) * cos(a_vy * PI * y / L) + 0.3e1 * a_uy * a_uy * u_y * u_y * pow(sin(a_uy * PI * y / L), 0.2e1) + 0.6e1 * a_uy * a_vx * u_y * v_x * sin(a_uy * PI * y / L) * sin(a_vx * PI * x / L) + 0.3e1 * a_vx * a_vx * v_x * v_x * pow(sin(a_vx * PI * x / L), 0.2e1) + 0.4e1 * a_vy * a_vy * v_y * v_y * pow(cos(a_vy * PI * y / L), 0.2e1)) * PI * PI * mu * pow(L, -0.2e1) / 0.3e1 - (0.4e1 * a_ux * a_ux * u_x * u_x * pow(cos(a_ux * PI * x / L), 0.2e1) - 0.4e1 * a_ux * a_vy * u_x * v_y * cos(a_ux * PI * x / L) * cos(a_vy * PI * y / L) + 0.3e1 * a_uy * a_uy * u_y * u_y * pow(sin(a_uy * PI * y / L), 0.2e1) + 0.6e1 * a_uy * a_vx * u_y * v_x * sin(a_uy * PI * y / L) * sin(a_vx * PI * x / L) + 0.3e1 * a_vx * a_vx * v_x * v_x * pow(sin(a_vx * PI * x / L), 0.2e1) + 0.4e1 * a_vy * a_vy * v_y * v_y * pow(cos(a_vy * PI * y / L), 0.2e1)) * PI * PI * mu_t * pow(L, -0.2e1) / 0.3e1 + cv * a_rhot * PI * rho_t * P * cos(a_rhot * PI * t / L) / L / R / RHO;
  return(Q_E);
}

// ----------------------------------------
//   Manufactured Analytical Solutions
// ----------------------------------------

template <typename Scalar>
Scalar MASA::fans_sa_transient_free_shear<Scalar>::eval_exact_u(Scalar x,Scalar y)
{
  Scalar u_an;
  u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  return u_an; 
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_free_shear<Scalar>::eval_exact_v(Scalar x,Scalar y)
{
  Scalar v_an;
  v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);
  return v_an;
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_free_shear<Scalar>::eval_exact_p(Scalar x,Scalar y)
{
  Scalar p_an;
  p_an = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L);
  return p_an;

}

template <typename Scalar>
Scalar MASA::fans_sa_transient_free_shear<Scalar>::eval_exact_nu(Scalar x,Scalar y)
{
  return eval_exact_nu(x,y,0.0);
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_free_shear<Scalar>::eval_exact_nu(Scalar x,Scalar y,Scalar t)
{
  Scalar nu_an;
  nu_an = nu_sa_0 + nu_sa_x * cos(a_nusax * pi * x / L) + nu_sa_y * cos(a_nusay * pi * y / L) + nu_sa_t * cos(a_nusat * pi * t / L);
  return nu_an;
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_free_shear<Scalar>::eval_exact_rho(Scalar x,Scalar y)
{
  Scalar rho_an;
  rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  return rho_an;
}




/* ------------------------------------------------
 *
 *   favre averaged navier stokes transient
 *            finite distance to wall
 *
 * -----------------------------------------------
 */ 


template <typename Scalar>
MASA::fans_sa_transient_d_finite<Scalar>::fans_sa_transient_d_finite()
{
  // change name to: fans_sa_transient_d_wall_bounded
  this->mmsname = "fans_sa_transient_d_finite";
  this->dimension=2;

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
  this->register_var("mu",&mu);

  this->register_var("u_t",&u_t);
  this->register_var("v_t",&v_t);
  this->register_var("p_t",&p_t);
  this->register_var("rho_t",&rho_t);
  this->register_var("a_ut",&a_ut);
  this->register_var("a_vt",&a_vt);
  this->register_var("a_pt",&a_pt);
  this->register_var("a_rhot",&a_rhot);

  this->register_var("nu_sa_0",&nu_sa_0);
  this->register_var("nu_sa_x",&nu_sa_x);
  this->register_var("nu_sa_y",&nu_sa_y);
  this->register_var("nu_sa_t",&nu_sa_y);

  this->register_var("a_nusax",&a_nusax);
  this->register_var("a_nusay",&a_nusay);
  this->register_var("a_nusat",&a_nusat);

  this->register_var("c_v1",&c_v1);
  this->register_var("c_b1",&c_b1);
  this->register_var("c_b2",&c_b2);
  this->register_var("c_w1",&c_w1);
  this->register_var("c_w2",&c_w2);
  this->register_var("c_w3",&c_w3);
  this->register_var("kappa",&kappa);
  this->register_var("sigma",&sigma);
  this->register_var("cp",&cp);
  this->register_var("cv",&cv);
  this->register_var("Pr",&Pr);
  this->register_var("Pr_t",&Pr_t);

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
  err += this->set_var("mu",.918);
  err += this->set_var("L",3.02);

  err += this->set_var("nu_sa_0",12.0);
  err += this->set_var("nu_sa_x",12.0);
  err += this->set_var("nu_sa_y",12.0);
  err += this->set_var("nu_sa_t",12.0);

  err += this->set_var("a_nusay",12.0);
  err += this->set_var("a_nusax",12.0);
  err += this->set_var("a_nusat",12.0);

  err += this->set_var("c_v1",12.0);
  err += this->set_var("c_b1",12.0);
  err += this->set_var("c_b2",12.0);
  err += this->set_var("c_w1",12.0);
  err += this->set_var("c_w2",12.0);
  err += this->set_var("c_w3",12.0);
  err += this->set_var("kappa",12.0);
  err += this->set_var("sigma",12.0);
  err += this->set_var("cp",12.0);
  err += this->set_var("cv",12.0);
  err += this->set_var("Pr",12.0);
  err += this->set_var("Pr_t",12.0);

  err += this->set_var("u_t",12.0);
  err += this->set_var("v_t",12.0);
  err += this->set_var("p_t",12.0);
  err += this->set_var("rho_t",12.0);
  err += this->set_var("a_ut",12.0);
  err += this->set_var("a_vt",12.0);
  err += this->set_var("a_pt",12.0);
  err += this->set_var("a_rhot",12.0);

  return err;
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_d_finite<Scalar>::eval_q_rho_u(Scalar x,Scalar y,Scalar t)
{
  Scalar Q_u;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar NU_SA;
  Scalar f_v1;
  Scalar chi;
  NU_SA = nu_sa_0 + nu_sa_x * cos(a_nusax * PI * x / L) + nu_sa_y * cos(a_nusay * PI * y / L) + nu_sa_t * cos(a_nusat * PI * t / L);
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  chi = RHO * NU_SA / mu;
  f_v1 = pow(chi, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1));
  Q_u = a_rhox * PI * rho_x * U * U * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * U * V * sin(a_rhoy * PI * y / L) / L - a_uy * PI * u_y * RHO * V * sin(a_uy * PI * y / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * U / L - a_px * PI * p_x * sin(a_px * PI * x / L) / L + (0.4e1 / 0.3e1 * a_ux * a_ux * u_x * sin(a_ux * PI * x / L) + a_uy * a_uy * u_y * cos(a_uy * PI * y / L)) * f_v1 * PI * PI * RHO * NU_SA * pow(L, -0.2e1) + (0.4e1 / 0.3e1 * a_ux * a_nusax * u_x * nu_sa_x * cos(a_ux * PI * x / L) * sin(a_nusax * PI * x / L) - a_uy * a_nusay * u_y * nu_sa_y * sin(a_uy * PI * y / L) * sin(a_nusay * PI * y / L) - a_vx * a_nusay * v_x * nu_sa_y * sin(a_vx * PI * x / L) * sin(a_nusay * PI * y / L) - 0.2e1 / 0.3e1 * a_vy * a_nusax * v_y * nu_sa_x * cos(a_vy * PI * y / L) * sin(a_nusax * PI * x / L)) * f_v1 * PI * PI * RHO * pow(L, -0.2e1) + (-0.4e1 / 0.3e1 * a_rhox * a_ux * rho_x * u_x * cos(a_rhox * PI * x / L) * cos(a_ux * PI * x / L) + 0.2e1 / 0.3e1 * a_rhox * a_vy * rho_x * v_y * cos(a_rhox * PI * x / L) * cos(a_vy * PI * y / L) - a_rhoy * a_uy * rho_y * u_y * sin(a_rhoy * PI * y / L) * sin(a_uy * PI * y / L) - a_rhoy * a_vx * rho_y * v_x * sin(a_rhoy * PI * y / L) * sin(a_vx * PI * x / L)) * f_v1 * PI * PI * NU_SA * pow(L, -0.2e1) + (pow(c_v1, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1)) + f_v1) * (0.4e1 * a_ux * a_ux * u_x * sin(a_ux * PI * x / L) + 0.3e1 * a_uy * a_uy * u_y * cos(a_uy * PI * y / L)) * PI * PI * mu * pow(L, -0.2e1) / 0.3e1 + a_rhot * PI * rho_t * U * cos(a_rhot * PI * t / L) / L - a_ut * PI * u_t * RHO * sin(a_ut * PI * t / L) / L;
  return(Q_u);

}

template <typename Scalar>
Scalar MASA::fans_sa_transient_d_finite<Scalar>::eval_q_rho_v(Scalar x,Scalar y,Scalar t)
{
  Scalar Q_v;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar NU_SA;
  Scalar f_v1;
  Scalar chi;
  NU_SA = nu_sa_0 + nu_sa_x * cos(a_nusax * PI * x / L) + nu_sa_y * cos(a_nusay * PI * y / L) + nu_sa_t * cos(a_nusat * PI * t / L);
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  chi = RHO * NU_SA / mu;
  f_v1 = pow(chi, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1));
  Q_v = a_rhox * PI * rho_x * U * V * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * V * sin(a_rhoy * PI * y / L) / L - a_vx * PI * v_x * RHO * U * sin(a_vx * PI * x / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * V / L + a_py * PI * p_y * cos(a_py * PI * y / L) / L + (a_vx * a_vx * v_x * cos(a_vx * PI * x / L) + 0.4e1 / 0.3e1 * a_vy * a_vy * v_y * sin(a_vy * PI * y / L)) * f_v1 * PI * PI * RHO * NU_SA * pow(L, -0.2e1) + (-0.2e1 / 0.3e1 * a_ux * a_nusay * u_x * nu_sa_y * cos(a_ux * PI * x / L) * sin(a_nusay * PI * y / L) - a_uy * a_nusax * u_y * nu_sa_x * sin(a_uy * PI * y / L) * sin(a_nusax * PI * x / L) - a_vx * a_nusax * v_x * nu_sa_x * sin(a_vx * PI * x / L) * sin(a_nusax * PI * x / L) + 0.4e1 / 0.3e1 * a_vy * a_nusay * v_y * nu_sa_y * cos(a_vy * PI * y / L) * sin(a_nusay * PI * y / L)) * f_v1 * PI * PI * RHO * pow(L, -0.2e1) + (a_rhox * a_uy * rho_x * u_y * cos(a_rhox * PI * x / L) * sin(a_uy * PI * y / L) + a_rhox * a_vx * rho_x * v_x * cos(a_rhox * PI * x / L) * sin(a_vx * PI * x / L) - 0.2e1 / 0.3e1 * a_rhoy * a_ux * rho_y * u_x * sin(a_rhoy * PI * y / L) * cos(a_ux * PI * x / L) + 0.4e1 / 0.3e1 * a_rhoy * a_vy * rho_y * v_y * sin(a_rhoy * PI * y / L) * cos(a_vy * PI * y / L)) * f_v1 * PI * PI * NU_SA * pow(L, -0.2e1) + (pow(c_v1, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1)) + f_v1) * (0.3e1 * a_vx * a_vx * v_x * cos(a_vx * PI * x / L) + 0.4e1 * a_vy * a_vy * v_y * sin(a_vy * PI * y / L)) * PI * PI * mu * pow(L, -0.2e1) / 0.3e1 + a_rhot * PI * rho_t * V * cos(a_rhot * PI * t / L) / L + a_vt * PI * v_t * RHO * cos(a_vt * PI * t / L) / L;
  return(Q_v);

}

template <typename Scalar>
Scalar MASA::fans_sa_transient_d_finite<Scalar>::eval_q_rho(Scalar x,Scalar y,Scalar t)
{
  Scalar Q_rho;
  Scalar RHO;
  Scalar U;
  Scalar V;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Q_rho = a_rhox * PI * rho_x * U * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * sin(a_rhoy * PI * y / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO / L + a_rhot * PI * rho_t * cos(a_rhot * PI * t / L) / L;
  return(Q_rho);

}

template <typename Scalar>
Scalar MASA::fans_sa_transient_d_finite<Scalar>::eval_q_nu(Scalar x,Scalar y,Scalar t)
{
  Scalar Q_nu;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar NU_SA;
  Scalar chi;
  Scalar f_v1;
  Scalar f_v2;
  Scalar Omega;
  Scalar Ssa;
  Scalar r;
  Scalar g;
  Scalar f_w;
  NU_SA = nu_sa_0 + nu_sa_x * cos(a_nusax * PI * x / L) + nu_sa_y * cos(a_nusay * PI * y / L) + nu_sa_t * cos(a_nusat * PI * t / L);
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  chi = NU_SA * RHO / mu;
  f_v1 = pow(chi, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1));
  f_v2 = 0.1e1 - chi / (0.1e1 + chi * f_v1);
  Omega = PI * sqrt(pow(a_uy * u_y * sin(a_uy * PI * y / L) - a_vx * v_x * sin(a_vx * PI * x / L), 0.2e1) * pow(L, -0.2e1));
  Ssa = Omega + NU_SA * f_v2 * pow(kappa, -0.2e1) * pow(y, -0.2e1);
  r = NU_SA / Ssa * pow(kappa, -0.2e1) * pow(y, -0.2e1);
  g = r + c_w2 * (pow(r, 0.6e1) - r);
  f_w = g * pow((0.1e1 + pow(c_w3, 0.6e1)) / (pow(g, 0.6e1) + pow(c_w3, 0.6e1)), 0.1e1 / 0.6e1);
  Q_nu = (a_nusax * a_nusax * nu_sa_x * cos(a_nusax * PI * x / L) + a_nusay * a_nusay * nu_sa_y * cos(a_nusay * PI * y / L)) * PI * PI * (NU_SA * RHO + mu) / sigma * pow(L, -0.2e1) - (a_nusax * a_nusax * nu_sa_x * nu_sa_x * pow(sin(a_nusax * PI * x / L), 0.2e1) + a_nusay * a_nusay * nu_sa_y * nu_sa_y * pow(sin(a_nusay * PI * y / L), 0.2e1)) * PI * PI * RHO * (c_b2 + 0.1e1) / sigma * pow(L, -0.2e1) + a_rhox * PI * rho_x * U * NU_SA * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * NU_SA * sin(a_rhoy * PI * y / L) / L - a_nusax * PI * nu_sa_x * RHO * U * sin(a_nusax * PI * x / L) / L - a_nusay * PI * nu_sa_y * RHO * V * sin(a_nusay * PI * y / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * NU_SA / L + (a_rhox * a_nusax * rho_x * nu_sa_x * cos(a_rhox * PI * x / L) * sin(a_nusax * PI * x / L) - a_rhoy * a_nusay * rho_y * nu_sa_y * sin(a_rhoy * PI * y / L) * sin(a_nusay * PI * y / L)) * PI * PI * NU_SA / sigma * pow(L, -0.2e1) + a_rhot * PI * rho_t * NU_SA * cos(a_rhot * PI * t / L) / L - a_nusat * PI * nu_sa_t * RHO * sin(a_nusat * PI * t / L) / L - c_b1 * Ssa * RHO * NU_SA + c_w1 * RHO * NU_SA * NU_SA * f_w * pow(y, -0.2e1);
  return(Q_nu);
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_d_finite<Scalar>::eval_q_rho_e(Scalar x,Scalar y,Scalar t)
{
  Scalar Q_E;
  /*

  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar P;
  Scalar NU_SA;
  Scalar chi;
  Scalar f_v1;
  Scalar f_v2;
  Scalar R;
  Scalar r;
  Scalar mu_t;
  Scalar Ssa;
  Scalar Omega;

  NU_SA = nu_sa_0 + nu_sa_x * cos(a_nusax * PI * x / L) + nu_sa_y * cos(a_nusay * PI * y / L) + nu_sa_t * cos(a_nusat * PI * t / L);
  r = NU_SA / Ssa * pow(kappa, -0.2e1) * pow(y, -0.2e1);
  RHO = r * rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / L);
  chi = RHO * NU_SA / mu;
  Ssa = Omega + NU_SA * f_v2 * pow(kappa, -0.2e1) * pow(y, -0.2e1);
  f_v1 = pow(chi, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1));
  f_v2 = 0.1e1 - chi / (0.1e1 + chi * f_v1);
  Omega = PI * sqrt(pow(a_uy * u_y * sin(a_uy * PI * y / L) - a_vx * v_x * sin(a_vx * PI * x / L), 0.2e1) * pow(L, -0.2e1));
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  P = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_t * cos(a_pt * PI * t / L);
  mu_t = RHO * NU_SA * f_v1;
  R = cp - cv;

  Q_E = -(0.4e1 * a_ux * a_ux * u_x * u_x * pow(cos(a_ux * PI * x / L), 0.2e1) - 0.4e1 * a_ux * a_vy * u_x * v_y * cos(a_ux * PI * x / L) * cos(a_vy * PI * y / L) + 0.3e1 * a_uy * a_uy * u_y * u_y * pow(sin(a_uy * PI * y / L), 0.2e1) + 0.6e1 * a_uy * a_vx * u_y * v_x * sin(a_uy * PI * y / L) * sin(a_vx * PI * x / L) + 0.3e1 * a_vx * a_vx * v_x * v_x * pow(sin(a_vx * PI * x / L), 0.2e1) + 0.4e1 * a_vy * a_vy * v_y * v_y * pow(cos(a_vy * PI * y / L), 0.2e1)) * PI * PI * mu_t * pow(L, -0.2e1) / 0.3e1 - (mu_t / Pr_t + mu / Pr) * (-(a_px * a_px * p_x * cos(a_px * PI * x / L) + a_py * a_py * p_y * sin(a_py * PI * y / L)) * cp * PI * PI * pow(L, -0.2e1) / R / RHO + (a_rhox * a_rhox * rho_x * sin(a_rhox * PI * x / L) + a_rhoy * a_rhoy * rho_y * cos(a_rhoy * PI * y / L)) * cp * PI * PI * P * pow(L, -0.2e1) / R * pow(RHO, -0.2e1)) + (0.4e1 * a_ux * a_ux * u_x * sin(a_ux * PI * x / L) + 0.3e1 * a_uy * a_uy * u_y * cos(a_uy * PI * y / L)) * PI * PI * mu_t * U * pow(L, -0.2e1) / 0.3e1 + (0.3e1 * a_vx * a_vx * v_x * cos(a_vx * PI * x / L) + 0.4e1 * a_vy * a_vy * v_y * sin(a_vy * PI * y / L)) * PI * PI * mu_t * V * pow(L, -0.2e1) / 0.3e1 + (0.3e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * U * U / L / 0.2e1 + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.3e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * V * V / L / 0.2e1 - (f_v1 + pow(c_v1, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1))) * (0.4e1 * a_ux * a_ux * u_x * u_x * pow(cos(a_ux * PI * x / L), 0.2e1) - 0.4e1 * a_ux * a_vy * u_x * v_y * cos(a_ux * PI * x / L) * cos(a_vy * PI * y / L) + 0.3e1 * a_uy * a_uy * u_y * u_y * pow(sin(a_uy * PI * y / L), 0.2e1) + 0.6e1 * a_uy * a_vx * u_y * v_x * sin(a_uy * PI * y / L) * sin(a_vx * PI * x / L) + 0.3e1 * a_vx * a_vx * v_x * v_x * pow(sin(a_vx * PI * x / L), 0.2e1) + 0.4e1 * a_vy * a_vy * v_y * v_y * pow(cos(a_vy * PI * y / L), 0.2e1)) * PI * PI * mu * pow(L, -0.2e1) / 0.3e1 - a_ut * PI * u_t * RHO * U * sin(a_ut * PI * t / L) / L + (U * U + V * V) * a_rhox * PI * rho_x * U * cos(a_rhox * PI * x / L) / L / 0.2e1 + a_vt * PI * v_t * RHO * V * cos(a_vt * PI * t / L) / L - (U * U + V * V) * a_rhoy * PI * rho_y * V * sin(a_rhoy * PI * y / L) / L / 0.2e1 - (a_px * a_nusax * p_x * nu_sa_x * sin(a_px * PI * x / L) * sin(a_nusax * PI * x / L) - a_py * a_nusay * p_y * nu_sa_y * cos(a_py * PI * y / L) * sin(a_nusay * PI * y / L)) * cp * PI * PI * mu_t * pow(L, -0.2e1) / Pr_t / R / RHO / NU_SA - (a_rhox * a_nusax * rho_x * nu_sa_x * cos(a_rhox * PI * x / L) * sin(a_nusax * PI * x / L) - a_rhoy * a_nusay * rho_y * nu_sa_y * sin(a_rhoy * PI * y / L) * sin(a_nusay * PI * y / L)) * cp * PI * PI * mu_t * P * pow(L, -0.2e1) / Pr_t / R * pow(RHO, -0.2e1) / NU_SA - (a_rhox * a_px * rho_x * p_x * cos(a_rhox * PI * x / L) * sin(a_px * PI * x / L) + a_rhoy * a_py * rho_y * p_y * sin(a_rhoy * PI * y / L) * cos(a_py * PI * y / L)) * (Pr * mu_t + 0.2e1 * Pr_t * mu) * cp * PI * PI / Pr / Pr_t * pow(L, -0.2e1) / R * pow(RHO, -0.2e1) - (a_rhox * a_rhox * rho_x * rho_x * pow(cos(a_rhox * PI * x / L), 0.2e1) + a_rhoy * a_rhoy * rho_y * rho_y * pow(sin(a_rhoy * PI * y / L), 0.2e1)) * (Pr * mu_t + 0.2e1 * Pr_t * mu) * cp * PI * PI * P / Pr / Pr_t * pow(L, -0.2e1) / R * pow(RHO, -0.3e1) + cv * a_rhot * PI * rho_t * P * cos(a_rhot * PI * t / L) / L / R / RHO + cp * a_py * PI * p_y * V * cos(a_py * PI * y / L) / L / R - cp * a_px * PI * p_x * U * sin(a_px * PI * x / L) / L / R + (0.4e1 / 0.3e1 * a_ux * a_nusax * u_x * nu_sa_x * cos(a_ux * PI * x / L) * sin(a_nusax * PI * x / L) - a_uy * a_nusay * u_y * nu_sa_y * sin(a_uy * PI * y / L) * sin(a_nusay * PI * y / L) - a_vx * a_nusay * v_x * nu_sa_y * sin(a_vx * PI * x / L) * sin(a_nusay * PI * y / L) - 0.2e1 / 0.3e1 * a_vy * a_nusax * v_y * nu_sa_x * cos(a_vy * PI * y / L) * sin(a_nusax * PI * x / L)) * PI * PI * f_v1 * RHO * U * pow(L, -0.2e1) + (-0.2e1 / 0.3e1 * a_ux * a_nusay * u_x * nu_sa_y * cos(a_ux * PI * x / L) * sin(a_nusay * PI * y / L) - a_uy * a_nusax * u_y * nu_sa_x * sin(a_uy * PI * y / L) * sin(a_nusax * PI * x / L) - a_vx * a_nusax * v_x * nu_sa_x * sin(a_vx * PI * x / L) * sin(a_nusax * PI * x / L) + 0.4e1 / 0.3e1 * a_vy * a_nusay * v_y * nu_sa_y * cos(a_vy * PI * y / L) * sin(a_nusay * PI * y / L)) * PI * PI * f_v1 * RHO * V * pow(L, -0.2e1) + (-0.4e1 / 0.3e1 * a_rhox * a_ux * rho_x * u_x * cos(a_rhox * PI * x / L) * cos(a_ux * PI * x / L) + 0.2e1 / 0.3e1 * a_rhox * a_vy * rho_x * v_y * cos(a_rhox * PI * x / L) * cos(a_vy * PI * y / L) - a_rhoy * a_uy * rho_y * u_y * sin(a_rhoy * PI * y / L) * sin(a_uy * PI * y / L) - a_rhoy * a_vx * rho_y * v_x * sin(a_rhoy * PI * y / L) * sin(a_vx * PI * x / L)) * PI * PI * f_v1 * U * NU_SA * pow(L, -0.2e1) + (a_rhox * a_uy * rho_x * u_y * cos(a_rhox * PI * x / L) * sin(a_uy * PI * y / L) + a_rhox * a_vx * rho_x * v_x * cos(a_rhox * PI * x / L) * sin(a_vx * PI * x / L) - 0.2e1 / 0.3e1 * a_rhoy * a_ux * rho_y * u_x * sin(a_rhoy * PI * y / L) * cos(a_ux * PI * x / L) + 0.4e1 / 0.3e1 * a_rhoy * a_vy * rho_y * v_y * sin(a_rhoy * PI * y / L) * cos(a_vy * PI * y / L)) * PI * PI * f_v1 * V * NU_SA * pow(L, -0.2e1) + (U * U + V * V) * a_rhot * PI * rho_t * cos(a_rhot * PI * t / L) / L / 0.2e1 + (f_v1 + pow(c_v1, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1))) * (0.4e1 * a_ux * a_ux * u_x * sin(a_ux * PI * x / L) + 0.3e1 * a_uy * a_uy * u_y * cos(a_uy * PI * y / L)) * PI * PI * mu * U * pow(L, -0.2e1) / 0.3e1 + (f_v1 + pow(c_v1, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1))) * (0.3e1 * a_vx * a_vx * v_x * cos(a_vx * PI * x / L) + 0.4e1 * a_vy * a_vy * v_y * sin(a_vy * PI * y / L)) * PI * PI * mu * V * pow(L, -0.2e1) / 0.3e1 - (a_uy * u_y * sin(a_uy * PI * y / L) + a_vx * v_x * sin(a_vx * PI * x / L)) * PI * RHO * U * V / L + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * cp * PI * P / L / R;
  */
  Q_E=42.0*x*y*t;
  return(Q_E);
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
MASA_INSTANTIATE_ALL(MASA::fans_sa_transient_free_shear);
