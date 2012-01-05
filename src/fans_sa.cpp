// -*-c++-*-
//
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
using namespace std;

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
  NU_SA = nu_sa_0 + nu_sa_x * std::cos(a_nusax * PI * x / L) + nu_sa_y * std::cos(a_nusay * PI * y / L) + nu_sa_t * std::cos(a_nusat * PI * t / L);
  RHO = rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_t * std::sin(a_rhot * PI * t / L);
  U = u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_t * std::cos(a_ut * PI * t / L);
  V = v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_t * std::sin(a_vt * PI * t / L);
  chi = RHO * NU_SA / mu;
  f_v1 = pow(chi, Scalar(0.3e1)) / (pow(chi, Scalar(0.3e1)) + pow(c_v1, Scalar(0.3e1)));
  Q_u = a_rhox * PI * rho_x * U * U * std::cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * U * V * std::sin(a_rhoy * PI * y / L) / L - a_uy * PI * u_y * RHO * V * std::sin(a_uy * PI * y / L) / L + (Scalar(0.2e1) * a_ux * u_x * std::cos(a_ux * PI * x / L) + a_vy * v_y * std::cos(a_vy * PI * y / L)) * PI * RHO * U / L - a_px * PI * p_x * std::sin(a_px * PI * x / L) / L + (Scalar(0.4e1) / Scalar(0.3e1) * a_ux * a_ux * u_x * std::sin(a_ux * PI * x / L) + a_uy * a_uy * u_y * std::cos(a_uy * PI * y / L)) * f_v1 * PI * PI * RHO * NU_SA * pow(L, Scalar(-0.2e1)) + (Scalar(0.4e1) / Scalar(0.3e1) * a_ux * a_nusax * u_x * nu_sa_x * std::cos(a_ux * PI * x / L) * std::sin(a_nusax * PI * x / L) - a_uy * a_nusay * u_y * nu_sa_y * std::sin(a_uy * PI * y / L) * std::sin(a_nusay * PI * y / L) - a_vx * a_nusay * v_x * nu_sa_y * std::sin(a_vx * PI * x / L) * std::sin(a_nusay * PI * y / L) - Scalar(0.2e1) / Scalar(0.3e1) * a_vy * a_nusax * v_y * nu_sa_x * std::cos(a_vy * PI * y / L) * std::sin(a_nusax * PI * x / L)) * f_v1 * PI * PI * RHO * pow(L, Scalar(-0.2e1)) + (Scalar(-0.4e1) / Scalar(0.3e1) * a_rhox * a_ux * rho_x * u_x * std::cos(a_rhox * PI * x / L) * std::cos(a_ux * PI * x / L) + Scalar(0.2e1) / Scalar(0.3e1) * a_rhox * a_vy * rho_x * v_y * std::cos(a_rhox * PI * x / L) * std::cos(a_vy * PI * y / L) - a_rhoy * a_uy * rho_y * u_y * std::sin(a_rhoy * PI * y / L) * std::sin(a_uy * PI * y / L) - a_rhoy * a_vx * rho_y * v_x * std::sin(a_rhoy * PI * y / L) * std::sin(a_vx * PI * x / L)) * f_v1 * PI * PI * NU_SA * pow(L, Scalar(-0.2e1)) + (pow(c_v1, Scalar(0.3e1)) / (pow(chi, Scalar(0.3e1)) + pow(c_v1, Scalar(0.3e1))) + f_v1) * (Scalar(0.4e1) * a_ux * a_ux * u_x * std::sin(a_ux * PI * x / L) + Scalar(0.3e1) * a_uy * a_uy * u_y * std::cos(a_uy * PI * y / L)) * PI * PI * mu * pow(L, Scalar(-0.2e1)) / Scalar(0.3e1) + a_rhot * PI * rho_t * U * std::cos(a_rhot * PI * t / L) / L - a_ut * PI * u_t * RHO * std::sin(a_ut * PI * t / L) / L;
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
  NU_SA = nu_sa_0 + nu_sa_x * std::cos(a_nusax * PI * x / L) + nu_sa_y * std::cos(a_nusay * PI * y / L) + nu_sa_t * std::cos(a_nusat * PI * t / L);
  RHO = rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_t * std::sin(a_rhot * PI * t / L);
  U = u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_t * std::cos(a_ut * PI * t / L);
  V = v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_t * std::sin(a_vt * PI * t / L);
  chi = RHO * NU_SA / mu;
  f_v1 = pow(chi, Scalar(0.3e1)) / (pow(chi, Scalar(0.3e1)) + pow(c_v1, Scalar(0.3e1)));
  Q_v = a_rhox * PI * rho_x * U * V * std::cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * V * std::sin(a_rhoy * PI * y / L) / L - a_vx * PI * v_x * RHO * U * std::sin(a_vx * PI * x / L) / L + (a_ux * u_x * std::cos(a_ux * PI * x / L) + Scalar(0.2e1) * a_vy * v_y * std::cos(a_vy * PI * y / L)) * PI * RHO * V / L + a_py * PI * p_y * std::cos(a_py * PI * y / L) / L + (a_vx * a_vx * v_x * std::cos(a_vx * PI * x / L) + Scalar(0.4e1) / Scalar(0.3e1) * a_vy * a_vy * v_y * std::sin(a_vy * PI * y / L)) * f_v1 * PI * PI * RHO * NU_SA * pow(L, Scalar(-0.2e1)) + (Scalar(-0.2e1) / Scalar(0.3e1) * a_ux * a_nusay * u_x * nu_sa_y * std::cos(a_ux * PI * x / L) * std::sin(a_nusay * PI * y / L) - a_uy * a_nusax * u_y * nu_sa_x * std::sin(a_uy * PI * y / L) * std::sin(a_nusax * PI * x / L) - a_vx * a_nusax * v_x * nu_sa_x * std::sin(a_vx * PI * x / L) * std::sin(a_nusax * PI * x / L) + Scalar(0.4e1) / Scalar(0.3e1) * a_vy * a_nusay * v_y * nu_sa_y * std::cos(a_vy * PI * y / L) * std::sin(a_nusay * PI * y / L)) * f_v1 * PI * PI * RHO * pow(L, Scalar(-0.2e1)) + (a_rhox * a_uy * rho_x * u_y * std::cos(a_rhox * PI * x / L) * std::sin(a_uy * PI * y / L) + a_rhox * a_vx * rho_x * v_x * std::cos(a_rhox * PI * x / L) * std::sin(a_vx * PI * x / L) - Scalar(0.2e1) / Scalar(0.3e1) * a_rhoy * a_ux * rho_y * u_x * std::sin(a_rhoy * PI * y / L) * std::cos(a_ux * PI * x / L) + Scalar(0.4e1) / Scalar(0.3e1) * a_rhoy * a_vy * rho_y * v_y * std::sin(a_rhoy * PI * y / L) * std::cos(a_vy * PI * y / L)) * f_v1 * PI * PI * NU_SA * pow(L, Scalar(-0.2e1)) + (pow(c_v1, Scalar(0.3e1)) / (pow(chi, Scalar(0.3e1)) + pow(c_v1, Scalar(0.3e1))) + f_v1) * (Scalar(0.3e1) * a_vx * a_vx * v_x * std::cos(a_vx * PI * x / L) + Scalar(0.4e1) * a_vy * a_vy * v_y * std::sin(a_vy * PI * y / L)) * PI * PI * mu * pow(L, Scalar(-0.2e1)) / Scalar(0.3e1) + a_rhot * PI * rho_t * V * std::cos(a_rhot * PI * t / L) / L + a_vt * PI * v_t * RHO * std::cos(a_vt * PI * t / L) / L;
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
  RHO = rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_t * std::sin(a_rhot * PI * t / L);
  U = u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_t * std::cos(a_ut * PI * t / L);
  V = v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_t * std::sin(a_vt * PI * t / L);
  Q_rho = a_rhox * PI * rho_x * U * std::cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * std::sin(a_rhoy * PI * y / L) / L + (a_ux * u_x * std::cos(a_ux * PI * x / L) + a_vy * v_y * std::cos(a_vy * PI * y / L)) * PI * RHO / L + a_rhot * PI * rho_t * std::cos(a_rhot * PI * t / L) / L;
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
  NU_SA = nu_sa_0 + nu_sa_x * std::cos(a_nusax * PI * x / L) + nu_sa_y * std::cos(a_nusay * PI * y / L) + nu_sa_t * std::cos(a_nusat * PI * t / L);
  RHO = rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_t * std::sin(a_rhot * PI * t / L);
  U = u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_t * std::cos(a_ut * PI * t / L);
  V = v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_t * std::sin(a_vt * PI * t / L);
  Q_nu = a_rhox * PI * rho_x * U * NU_SA * std::cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * NU_SA * std::sin(a_rhoy * PI * y / L) / L - a_nusax * PI * nu_sa_x * RHO * U * std::sin(a_nusax * PI * x / L) / L - a_nusay * PI * nu_sa_y * RHO * V * std::sin(a_nusay * PI * y / L) / L + a_rhot * PI * rho_t * NU_SA * std::cos(a_rhot * PI * t / L) / L - a_nusat * PI * nu_sa_t * RHO * std::sin(a_nusat * PI * t / L) / L - (a_nusax * a_nusax * nu_sa_x * nu_sa_x * pow(std::sin(a_nusax * PI * x / L), Scalar(0.2e1)) + a_nusay * a_nusay * nu_sa_y * nu_sa_y * pow(std::sin(a_nusay * PI * y / L), Scalar(0.2e1))) * (Scalar(0.1e1) + c_b2) * PI * PI * RHO / sigma * pow(L, Scalar(-0.2e1)) + (a_rhox * a_nusax * rho_x * nu_sa_x * std::cos(a_rhox * PI * x / L) * std::sin(a_nusax * PI * x / L) - a_rhoy * a_nusay * rho_y * nu_sa_y * std::sin(a_rhoy * PI * y / L) * std::sin(a_nusay * PI * y / L)) * PI * PI * NU_SA / sigma * pow(L, Scalar(-0.2e1)) - c_b1 * std::sqrt(pow(-a_uy * u_y * std::sin(a_uy * PI * y / L) + a_vx * v_x * std::sin(a_vx * PI * x / L), Scalar(0.2e1)) * pow(L, Scalar(-0.2e1))) * PI * RHO * NU_SA + (a_ux * u_x * std::cos(a_ux * PI * x / L) + a_vy * v_y * std::cos(a_vy * PI * y / L)) * PI * RHO * NU_SA / L + (a_nusax * a_nusax * nu_sa_x * std::cos(a_nusax * PI * x / L) + a_nusay * a_nusay * nu_sa_y * std::cos(a_nusay * PI * y / L)) * (RHO * NU_SA + mu) * PI * PI / sigma * pow(L, Scalar(-0.2e1));
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
  Scalar Q_E;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar P;
  Scalar NU_SA;
  Scalar chi;
  Scalar f_v1;
  Scalar cp;
  Scalar cv;
  Scalar mu_t;
  
  NU_SA = nu_sa_0 + nu_sa_x * std::cos(a_nusax * PI * x / L) + nu_sa_y * std::cos(a_nusay * PI * y / L) + nu_sa_t * std::cos(a_nusat * PI * t / L);
  RHO = rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_t * std::sin(a_rhot * PI * t / L);
  U = u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_t * std::cos(a_ut * PI * t / L);
  V = v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_t * std::sin(a_vt * PI * t / L);
  P = p_0 + p_x * std::cos(a_px * PI * x / L) + p_y * std::sin(a_py * PI * y / L) + p_t * std::cos(a_pt * PI * t / L);
  chi = RHO * NU_SA / mu;
  f_v1 = pow(chi, Scalar(0.3e1)) / (pow(chi, Scalar(0.3e1)) + pow(c_v1, Scalar(0.3e1)));
  mu_t = RHO * NU_SA * f_v1;
  cv = R/(Gamma-1.0);
  cp = Gamma*cv;
  Q_E = -(mu_t / Pr_t + mu / Pr) * (-(a_px * a_px * p_x * std::cos(a_px * PI * x / L) + a_py * a_py * p_y * std::sin(a_py * PI * y / L)) * cp * PI * PI * pow(L, Scalar(-0.2e1)) / R / RHO + (a_rhox * a_rhox * rho_x * std::sin(a_rhox * PI * x / L) + a_rhoy * a_rhoy * rho_y * std::cos(a_rhoy * PI * y / L)) * cp * PI * PI * P * pow(L, Scalar(-0.2e1)) / R * pow(RHO, Scalar(-0.2e1))) + (U * U + V * V) * a_rhox * PI * rho_x * U * std::cos(a_rhox * PI * x / L) / L / Scalar(0.2e1) - (a_rhox * a_rhox * rho_x * rho_x * pow(std::cos(a_rhox * PI * x / L), Scalar(0.2e1)) + a_rhoy * a_rhoy * rho_y * rho_y * pow(std::sin(a_rhoy * PI * y / L), Scalar(0.2e1))) * (Pr * mu_t + Scalar(0.2e1) * Pr_t * mu) * cp * PI * PI * P / Pr / Pr_t * pow(L, Scalar(-0.2e1)) / R * pow(RHO, Scalar(-0.3e1)) - (a_rhox * a_px * rho_x * p_x * std::cos(a_rhox * PI * x / L) * std::sin(a_px * PI * x / L) + a_rhoy * a_py * rho_y * p_y * std::sin(a_rhoy * PI * y / L) * std::cos(a_py * PI * y / L)) * (Pr * mu_t + Scalar(0.2e1) * Pr_t * mu) * cp * PI * PI / Pr / Pr_t * pow(L, Scalar(-0.2e1)) / R * pow(RHO, Scalar(-0.2e1)) - (a_rhox * a_nusax * rho_x * nu_sa_x * std::cos(a_rhox * PI * x / L) * std::sin(a_nusax * PI * x / L) - a_rhoy * a_nusay * rho_y * nu_sa_y * std::sin(a_rhoy * PI * y / L) * std::sin(a_nusay * PI * y / L)) * cp * PI * PI * mu_t * P * pow(L, Scalar(-0.2e1)) / Pr_t / R * pow(RHO, Scalar(-0.2e1)) / NU_SA - (a_px * a_nusax * p_x * nu_sa_x * std::sin(a_px * PI * x / L) * std::sin(a_nusax * PI * x / L) - a_py * a_nusay * p_y * nu_sa_y * std::cos(a_py * PI * y / L) * std::sin(a_nusay * PI * y / L)) * cp * PI * PI * mu_t * pow(L, Scalar(-0.2e1)) / Pr_t / R / RHO / NU_SA + (Scalar(0.4e1) / Scalar(0.3e1) * a_ux * a_nusax * u_x * nu_sa_x * std::cos(a_ux * PI * x / L) * std::sin(a_nusax * PI * x / L) - a_uy * a_nusay * u_y * nu_sa_y * std::sin(a_uy * PI * y / L) * std::sin(a_nusay * PI * y / L) - a_vx * a_nusay * v_x * nu_sa_y * std::sin(a_vx * PI * x / L) * std::sin(a_nusay * PI * y / L) - Scalar(0.2e1) / Scalar(0.3e1) * a_vy * a_nusax * v_y * nu_sa_x * std::cos(a_vy * PI * y / L) * std::sin(a_nusax * PI * x / L)) * PI * PI * f_v1 * RHO * U * pow(L, Scalar(-0.2e1)) + (Scalar(-0.2e1) / Scalar(0.3e1) * a_ux * a_nusay * u_x * nu_sa_y * std::cos(a_ux * PI * x / L) * std::sin(a_nusay * PI * y / L) - a_uy * a_nusax * u_y * nu_sa_x * std::sin(a_uy * PI * y / L) * std::sin(a_nusax * PI * x / L) - a_vx * a_nusax * v_x * nu_sa_x * std::sin(a_vx * PI * x / L) * std::sin(a_nusax * PI * x / L) + Scalar(0.4e1) / Scalar(0.3e1) * a_vy * a_nusay * v_y * nu_sa_y * std::cos(a_vy * PI * y / L) * std::sin(a_nusay * PI * y / L)) * PI * PI * f_v1 * RHO * V * pow(L, Scalar(-0.2e1)) + (Scalar(-0.4e1) / Scalar(0.3e1) * a_rhox * a_ux * rho_x * u_x * std::cos(a_rhox * PI * x / L) * std::cos(a_ux * PI * x / L) + Scalar(0.2e1) / Scalar(0.3e1) * a_rhox * a_vy * rho_x * v_y * std::cos(a_rhox * PI * x / L) * std::cos(a_vy * PI * y / L) - a_rhoy * a_uy * rho_y * u_y * std::sin(a_rhoy * PI * y / L) * std::sin(a_uy * PI * y / L) - a_rhoy * a_vx * rho_y * v_x * std::sin(a_rhoy * PI * y / L) * std::sin(a_vx * PI * x / L)) * PI * PI * f_v1 * U * NU_SA * pow(L, Scalar(-0.2e1)) + (a_rhox * a_uy * rho_x * u_y * std::cos(a_rhox * PI * x / L) * std::sin(a_uy * PI * y / L) + a_rhox * a_vx * rho_x * v_x * std::cos(a_rhox * PI * x / L) * std::sin(a_vx * PI * x / L) - Scalar(0.2e1) / Scalar(0.3e1) * a_rhoy * a_ux * rho_y * u_x * std::sin(a_rhoy * PI * y / L) * std::cos(a_ux * PI * x / L) + Scalar(0.4e1) / Scalar(0.3e1) * a_rhoy * a_vy * rho_y * v_y * std::sin(a_rhoy * PI * y / L) * std::cos(a_vy * PI * y / L)) * PI * PI * f_v1 * V * NU_SA * pow(L, Scalar(-0.2e1)) + (U * U + V * V) * a_rhot * PI * rho_t * std::cos(a_rhot * PI * t / L) / L / Scalar(0.2e1) + (f_v1 + pow(c_v1, Scalar(0.3e1)) / (pow(chi, Scalar(0.3e1)) + pow(c_v1, Scalar(0.3e1)))) * (Scalar(0.4e1) * a_ux * a_ux * u_x * std::sin(a_ux * PI * x / L) + Scalar(0.3e1) * a_uy * a_uy * u_y * std::cos(a_uy * PI * y / L)) * PI * PI * mu * U * pow(L, Scalar(-0.2e1)) / Scalar(0.3e1) + (f_v1 + pow(c_v1, Scalar(0.3e1)) / (pow(chi, Scalar(0.3e1)) + pow(c_v1, Scalar(0.3e1)))) * (Scalar(0.3e1) * a_vx * a_vx * v_x * std::cos(a_vx * PI * x / L) + Scalar(0.4e1) * a_vy * a_vy * v_y * std::sin(a_vy * PI * y / L)) * PI * PI * mu * V * pow(L, Scalar(-0.2e1)) / Scalar(0.3e1) - (a_uy * u_y * std::sin(a_uy * PI * y / L) + a_vx * v_x * std::sin(a_vx * PI * x / L)) * PI * RHO * U * V / L + (a_ux * u_x * std::cos(a_ux * PI * x / L) + a_vy * v_y * std::cos(a_vy * PI * y / L)) * cp * PI * P / L / R - a_ut * PI * u_t * RHO * U * std::sin(a_ut * PI * t / L) / L - (U * U + V * V) * a_rhoy * PI * rho_y * V * std::sin(a_rhoy * PI * y / L) / L / Scalar(0.2e1) + a_vt * PI * v_t * RHO * V * std::cos(a_vt * PI * t / L) / L + cp * a_py * PI * p_y * V * std::cos(a_py * PI * y / L) / L / R - cp * a_px * PI * p_x * U * std::sin(a_px * PI * x / L) / L / R + (Scalar(0.4e1) * a_ux * a_ux * u_x * std::sin(a_ux * PI * x / L) + Scalar(0.3e1) * a_uy * a_uy * u_y * std::cos(a_uy * PI * y / L)) * PI * PI * mu_t * U * pow(L, Scalar(-0.2e1)) / Scalar(0.3e1) + (Scalar(0.3e1) * a_vx * a_vx * v_x * std::cos(a_vx * PI * x / L) + Scalar(0.4e1) * a_vy * a_vy * v_y * std::sin(a_vy * PI * y / L)) * PI * PI * mu_t * V * pow(L, Scalar(-0.2e1)) / Scalar(0.3e1) + (Scalar(0.3e1) * a_ux * u_x * std::cos(a_ux * PI * x / L) + a_vy * v_y * std::cos(a_vy * PI * y / L)) * PI * RHO * U * U / L / Scalar(0.2e1) + (a_ux * u_x * std::cos(a_ux * PI * x / L) + Scalar(0.3e1) * a_vy * v_y * std::cos(a_vy * PI * y / L)) * PI * RHO * V * V / L / Scalar(0.2e1) - (f_v1 + pow(c_v1, Scalar(0.3e1)) / (pow(chi, Scalar(0.3e1)) + pow(c_v1, Scalar(0.3e1)))) * (Scalar(0.4e1) * a_ux * a_ux * u_x * u_x * pow(std::cos(a_ux * PI * x / L), Scalar(0.2e1)) - Scalar(0.4e1) * a_ux * a_vy * u_x * v_y * std::cos(a_ux * PI * x / L) * std::cos(a_vy * PI * y / L) + Scalar(0.3e1) * a_uy * a_uy * u_y * u_y * pow(std::sin(a_uy * PI * y / L), Scalar(0.2e1)) + Scalar(0.6e1) * a_uy * a_vx * u_y * v_x * std::sin(a_uy * PI * y / L) * std::sin(a_vx * PI * x / L) + Scalar(0.3e1) * a_vx * a_vx * v_x * v_x * pow(std::sin(a_vx * PI * x / L), Scalar(0.2e1)) + Scalar(0.4e1) * a_vy * a_vy * v_y * v_y * pow(std::cos(a_vy * PI * y / L), Scalar(0.2e1))) * PI * PI * mu * pow(L, Scalar(-0.2e1)) / Scalar(0.3e1) - (Scalar(0.4e1) * a_ux * a_ux * u_x * u_x * pow(std::cos(a_ux * PI * x / L), Scalar(0.2e1)) - Scalar(0.4e1) * a_ux * a_vy * u_x * v_y * std::cos(a_ux * PI * x / L) * std::cos(a_vy * PI * y / L) + Scalar(0.3e1) * a_uy * a_uy * u_y * u_y * pow(std::sin(a_uy * PI * y / L), Scalar(0.2e1)) + Scalar(0.6e1) * a_uy * a_vx * u_y * v_x * std::sin(a_uy * PI * y / L) * std::sin(a_vx * PI * x / L) + Scalar(0.3e1) * a_vx * a_vx * v_x * v_x * pow(std::sin(a_vx * PI * x / L), Scalar(0.2e1)) + Scalar(0.4e1) * a_vy * a_vy * v_y * v_y * pow(std::cos(a_vy * PI * y / L), Scalar(0.2e1))) * PI * PI * mu_t * pow(L, Scalar(-0.2e1)) / Scalar(0.3e1) + cv * a_rhot * PI * rho_t * P * std::cos(a_rhot * PI * t / L) / L / R / RHO;
  return(Q_E);
}

// ----------------------------------------
//   Manufactured Analytical Solutions
// ----------------------------------------

template <typename Scalar>
Scalar MASA::fans_sa_transient_free_shear<Scalar>::eval_exact_u(Scalar x,Scalar y)
{
  Scalar u_an;
  u_an = u_0 + u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L);
  return u_an; 
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_free_shear<Scalar>::eval_exact_v(Scalar x,Scalar y)
{
  Scalar v_an;
  v_an = v_0 + v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L);
  return v_an;
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_free_shear<Scalar>::eval_exact_p(Scalar x,Scalar y)
{
  Scalar p_an;
  p_an = p_0 + p_x * std::cos(a_px * pi * x / L) + p_y * std::sin(a_py * pi * y / L);
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
  nu_an = nu_sa_0 + nu_sa_x * std::cos(a_nusax * pi * x / L) + nu_sa_y * std::cos(a_nusay * pi * y / L) + nu_sa_t * std::cos(a_nusat * pi * t / L);
  return nu_an;
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_free_shear<Scalar>::eval_exact_rho(Scalar x,Scalar y)
{
  Scalar rho_an;
  rho_an = rho_0 + rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L);
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
MASA::fans_sa_steady_wall_bounded<Scalar>::fans_sa_steady_wall_bounded()
{
  this->mmsname = "fans_sa_steady_wall_bounded";
  this->dimension=2;

  // fluid properties
  this->register_var("mu",&mu);
  this->register_var("R",&R);
  this->register_var("p_0",&p_0);
  this->register_var("Pr",&Pr);
  this->register_var("Pr_t",&Pr_t);

  // SA calibration model
  this->register_var("eta1",&eta1);
  this->register_var("eta_v",&eta_v);
  this->register_var("kappa",&kappa);
  this->register_var("sigma",&sigma);
  this->register_var("c_b1",&c_b1);
  this->register_var("c_b2",&c_b2);
  this->register_var("c_v1",&c_v1);
  this->register_var("c_v2",&c_v2);
  this->register_var("c_v3",&c_v3);
  this->register_var("c_w2",&c_w2);
  this->register_var("c_w3",&c_w3);

  // manufactured solutions
  this->register_var("T_inf",&T_inf);
  this->register_var("M_inf",&M_inf);
  this->register_var("r_T",&r_T);
  this->register_var("Gamma",&Gamma);
  this->register_var("alpha",&alpha);
  this->register_var("C_cf",&C_cf);
  this->register_var("C",&C);
  this->register_var("b",&b);

  // initalize
  this->init_var();
  
}//done with constructor

template <typename Scalar>
int MASA::fans_sa_steady_wall_bounded<Scalar>::init_var()
{
  int err = 0;

  // fluid properties
  err += this->set_var("mu",0.1);//0.0001);
  err += this->set_var("R",287);
  err += this->set_var("p_0",1e3);//1e4); //this->set_var("p_0",1000);
  err += this->set_var("Pr",0.71);
  err += this->set_var("Pr_t",0.9);

  // SA calibration model
  err += this->set_var("eta1",11.0);
  err += this->set_var("eta_v",30.0);
  err += this->set_var("kappa",0.41);
  err += this->set_var("sigma",2.0/3.0);
  err += this->set_var("c_b1",0.1355);
  err += this->set_var("c_b2",0.622);
  err += this->set_var("c_v1",7.1);
  err += this->set_var("c_v2",0.7); //12);
  err += this->set_var("c_v3",0.9); //12);
  err += this->set_var("c_w2",0.3);
  err += this->set_var("c_w3",2.0);

  // manufactured solutions
  err += this->set_var("T_inf",250.0);
  err += this->set_var("M_inf",0.8);
  err += this->set_var("r_T",0.9);
  err += this->set_var("Gamma",1.4);
  err += this->set_var("alpha",5.0);
  err += this->set_var("C_cf",0.027);  
  err += this->set_var("C",5.0);
  err += this->set_var("b",0.33);

  // this sets all constants set from parameters
  update(0,0);

  return err;
}

template <typename Scalar>
Scalar MASA::fans_sa_steady_wall_bounded<Scalar>::eval_q_rho_u(Scalar x,Scalar y)
{
  // this resets all constants set from parameters
  update(x,y);

  // "Contribution from the convective terms to the total source term ---------------------------------------"
  Scalar Q_u_convection = r_T * (Gamma - 0.1e1) * M_inf * M_inf * T_inf * y_plus * u_tau * d_ueqplus_yplus * RHO * U * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / y / T - r_T * (Gamma - 0.1e1) * (u_eq_plus + y_plus * d_ueqplus_yplus) * M_inf * M_inf * T_inf * u_tau * RHO * pow(U, 0.3e1) * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / x / T / 0.14e2 + y_plus * u_tau * d_ueqplus_yplus * RHO * V * cos(A / u_inf * u_eq) / y - (u_eq_plus + y_plus * d_ueqplus_yplus) * u_tau * RHO * U * cos(A / u_inf * u_eq) / x / 0.7e1 + RHO * U * V / y;

  // "Contribution from the gradient of pressure (body forces) to the total source term -----------------"
  Scalar Q_u_gradp = 0.0e0;
  
  // "Contribution  from the viscous terms to the total source term -----------------------------------"
  Scalar Q_u_viscous = -(Gamma - 0.1e1) * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * r_T * mu_t * M_inf * M_inf * T_inf * y_plus * y_plus * u_tau * u_tau * d_ueqplus_yplus * d_ueqplus_yplus * U * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(u_inf, -0.2e1) * pow(y, -0.2e1) / NU_SA / RHO / T - pow(u_eq_plus + y_plus * d_ueqplus_yplus, 0.2e1) * (Gamma - 0.1e1) * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * r_T * mu_t * M_inf * M_inf * T_inf * u_tau * u_tau * U * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(x, -0.2e1) * pow(u_inf, -0.2e1) / NU_SA / RHO / T / 0.147e3 + (Gamma - 0.1e1) * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (0.43e2 * y_plus * d_ueqplus_yplus - 0.2e1 * u_eq_plus) * r_T * mu_t * M_inf * M_inf * T_inf * u_tau * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / x / y / NU_SA / RHO / T / 0.42e2 + (0.2e1 * alpha * y - kappa * u_tau) * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * mu_t * y_plus * u_tau * d_ueqplus_yplus * cos(A / u_inf * u_eq) / y * pow(NU_SA, -0.2e1) / RHO - (u_eq_plus + y_plus * d_ueqplus_yplus) * kappa * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * mu_t * y * u_tau * u_tau * cos(A / u_inf * u_eq) * pow(x, -0.2e1) * pow(NU_SA, -0.2e1) / RHO / 0.147e3 - (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (0.90e2 * alpha * y - 0.43e2 * kappa * u_tau) * mu_t * V / x * pow(NU_SA, -0.2e1) / RHO / 0.42e2 + (-0.2e1 * mu_t - 0.2e1 * mu) * (0.2e1 / 0.3e1 * D2uDx2 + D2vDxy / 0.6e1 + D2uDy2 / 0.2e1);

    // "Total source term --------------------------------------------------------------------------------"
  Scalar Q_u = Q_u_convection + Q_u_gradp + Q_u_viscous;
  return(Q_u);
}

template <typename Scalar>
Scalar MASA::fans_sa_steady_wall_bounded<Scalar>::eval_q_rho_v(Scalar x,Scalar y)
{
  // this resets all constants set from parameters
  update(x,y);

  // "Contribution from the convective terms to the total source term -------------------------------------------"
  Scalar Q_v_convection = r_T * (Gamma - 0.1e1) * M_inf * M_inf * T_inf * y_plus * u_tau * d_ueqplus_yplus * RHO * U * V * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / y / T - r_T * (Gamma - 0.1e1) * (y_plus * d_ueqplus_yplus + u_eq_plus) * M_inf * M_inf * T_inf * u_tau * RHO * U * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / x / T / 0.14e2 - (y_plus * d_ueqplus_yplus + u_eq_plus) * u_tau * RHO * V * cos(A / u_inf * u_eq) / x / 0.14e2 - 0.15e2 / 0.14e2 * RHO * U * V / x + 0.2e1 * RHO * V * V / y;

  // "Contribution from the gradient of pressure (body forces) to the total source term ------------------------"
  Scalar Q_v_gradp = 0.0e0;

  // "Contribution  from the viscous terms to the total source term ----------------------------------------------"
  Scalar Q_v_viscous = (Gamma - 0.1e1) * (y_plus * d_ueqplus_yplus + u_eq_plus) * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * r_T * mu_t * M_inf * M_inf * T_inf * u_tau * u_tau * y_plus * U * d_ueqplus_yplus * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(u_inf, -0.2e1) / x / y / NU_SA / RHO / T / 0.42e2 - 0.4e1 / 0.3e1 * (Gamma - 0.1e1) * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * r_T * mu_t * M_inf * M_inf * T_inf * u_tau * y_plus * U * V * d_ueqplus_yplus * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) * pow(y, -0.2e1) / NU_SA / RHO / T - 0.15e2 / 0.196e3 * (Gamma - 0.1e1) * (y_plus * d_ueqplus_yplus + u_eq_plus) * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * r_T * mu_t * M_inf * M_inf * T_inf * u_tau * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) * pow(x, -0.2e1) / NU_SA / RHO / T + kappa * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * mu_t * u_tau * u_tau * y_plus * d_ueqplus_yplus * cos(A / u_inf * u_eq) / x * pow(NU_SA, -0.2e1) / RHO / 0.14e2 + (y_plus * d_ueqplus_yplus + u_eq_plus) * (0.2e1 * alpha * y - kappa * u_tau) * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * mu_t * u_tau * cos(A / u_inf * u_eq) / x * pow(NU_SA, -0.2e1) / RHO / 0.21e2 - 0.15e2 / 0.196e3 * kappa * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * mu_t * y * u_tau * V * pow(x, -0.2e1) * pow(NU_SA, -0.2e1) / RHO + 0.4e1 / 0.3e1 * (0.2e1 * alpha * y - kappa * u_tau) * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * mu_t * V / y * pow(NU_SA, -0.2e1) / RHO + (-0.2e1 * mu - 0.2e1 * mu_t) * (D2uDxy / 0.6e1 + D2vDx2 / 0.2e1 + 0.2e1 / 0.3e1 * D2vDy2);

    // "Total source term ----------------------------------------------------------------------------------"
  Scalar Q_v = Q_v_convection + Q_v_gradp + Q_v_viscous;
  return(Q_v);
}

template <typename Scalar>
Scalar MASA::fans_sa_steady_wall_bounded<Scalar>::eval_q_rho(Scalar x,Scalar y)
{
  // this resets all constants set from parameters
  update(x,y);
  
  // "Contribution from the convective terms to the total source term ---------------------------"
  Scalar Q_rho_convection = r_T * (Gamma - 0.1e1) * d_ueqplus_yplus * M_inf * M_inf * T_inf * y_plus * u_tau * RHO * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / y / T - r_T * (Gamma - 0.1e1) * (d_ueqplus_yplus * y_plus + u_eq_plus) * M_inf * M_inf * T_inf * u_tau * RHO * U * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / x / T / 0.14e2 - (d_ueqplus_yplus * y_plus + u_eq_plus) * u_tau * RHO * cos(A / u_inf * u_eq) / x / 0.14e2 + RHO * V / y;

  // "Total source term -------------------------------------------------------------------------"
  Scalar Q_rho = Q_rho_convection;
  return(Q_rho);
}

template <typename Scalar>
Scalar MASA::fans_sa_steady_wall_bounded<Scalar>::eval_q_nu(Scalar x,Scalar y)
{
  // this resets all constants set from parameters
  update(x,y);

  // "Contribution from the convective terms to the total source term --------------------------------------"
  Scalar Q_nusa_convection = r_T * (Gamma - 0.1e1) * d_ueqplus_yplus * M_inf * M_inf * T_inf * y_plus * u_tau * NU_SA * RHO * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / y / T - r_T * (Gamma - 0.1e1) * (d_ueqplus_yplus * y_plus + u_eq_plus) * M_inf * M_inf * T_inf * u_tau * NU_SA * RHO * U * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / x / T / 0.14e2 - kappa * u_tau * y * RHO * U / x / 0.14e2 - (d_ueqplus_yplus * y_plus + u_eq_plus) * u_tau * NU_SA * RHO * cos(A / u_inf * u_eq) / x / 0.14e2 + RHO * NU_SA * V / y + (-0.2e1 * alpha * y + kappa * u_tau) * RHO * V;
  
  // "Contribution from theproduction  to the total source term ------------------------------------------------"
  Scalar Q_nusa_production = -c_b1 * (Sm + Omega) * RHO * NU_SA;

  // "Contribution  from the diffusion to the total source term ----------------------------------------------"
  Scalar Q_nusa_diffusion = (0.2e1 * alpha * y - kappa * u_tau) * r_T * (Gamma - 0.1e1) * M_inf * M_inf * T_inf * y_plus * u_tau * d_ueqplus_yplus * NU_SA * RHO * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / sigma / y / T - r_T * (Gamma - 0.1e1) * (d_ueqplus_yplus * y_plus + u_eq_plus) * kappa * y * M_inf * M_inf * T_inf * u_tau * u_tau * NU_SA * RHO * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / sigma * pow(x, -0.2e1) / T / 0.196e3 - kappa * kappa * y * y * u_tau * u_tau * RHO / sigma * pow(x, -0.2e1) / 0.196e3 - pow(0.2e1 * alpha * y - kappa * u_tau, 0.2e1) * RHO / sigma - 0.15e2 / 0.196e3 * kappa * (RHO * NU_SA + mu) * y * u_tau / sigma * pow(x, -0.2e1) + 0.2e1 * alpha * (RHO * NU_SA + mu) / sigma;

  // "Contribution  from the grad squared  to the total source term ----------------------------------------"
  Scalar Q_nusa_gradsquare = -c_b2 * kappa * kappa * y * y * u_tau * u_tau * RHO / sigma * pow(x, -0.2e1) / 0.196e3 - c_b2 * pow(0.2e1 * alpha * y - kappa * u_tau, 0.2e1) * RHO / sigma;

  // "Contribution  from the dissipation to the total source term ----------------------------------------------"
  Scalar Q_nusa_dissipation = c_w1 * f_w * RHO * NU_SA * NU_SA * pow(d, Scalar(-0.2e1));

  // "Total source term ----------------------------------------------------------------------------------"
  Scalar Q_nu_sa = Q_nusa_convection + Q_nusa_production + Q_nusa_diffusion + Q_nusa_gradsquare + Q_nusa_dissipation;
  return(Q_nu_sa);

}

template <typename Scalar>
Scalar MASA::fans_sa_steady_wall_bounded<Scalar>::eval_q_rho_e(Scalar x,Scalar y)
{
  // this resets all constants set from parameters
  update(x,y);

  // "Contribution from the convection to the total source term --------------------------------------------------"
  Scalar Q_E_convection = T_inf * r_T * (Gamma - 0.1e1) * M_inf * M_inf * y_plus * u_tau * d_ueqplus_yplus * RHO * pow(U, 0.3e1) * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / y / T / 0.2e1 - r_T * (Gamma - 0.1e1) * (u_eq_plus + y_plus * d_ueqplus_yplus) * (U * U + V * V) * M_inf * M_inf * T_inf * u_tau * RHO * U * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / x / T / 0.28e2 + (r_T * Gamma * M_inf * M_inf * T_inf * V * V - r_T * M_inf * M_inf * T_inf * V * V + 0.2e1 * u_inf * u_inf * T) * y_plus * u_tau * d_ueqplus_yplus * RHO * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / y / T / 0.2e1 - 0.15e2 / 0.14e2 * RHO * U * V * V / x - (u_eq_plus + y_plus * d_ueqplus_yplus) * (0.3e1 * U * U + V * V + 0.2e1 * cp * T) * u_tau * RHO * cos(A / u_inf * u_eq) / x / 0.28e2 + (U * U + 0.3e1 * V * V + 0.2e1 * cp * T) * RHO * V / y / 0.2e1;

  // "Contribution from the heat flux to the total source term --------------------------------------------------"
  Scalar Q_E_heat_flux = (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * pow(Gamma - 0.1e1, 0.2e1) * cp * r_T * r_T * mu_t * pow(M_inf, 0.4e1) * T_inf * T_inf * y_plus * y_plus * u_tau * u_tau * d_ueqplus_yplus * d_ueqplus_yplus * U * U * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(u_inf, -0.4e1) / Pr_t * pow(y, -0.2e1) / NU_SA / RHO / T + (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * pow(Gamma - 0.1e1, 0.2e1) * pow(u_eq_plus + y_plus * d_ueqplus_yplus, 0.2e1) * cp * r_T * r_T * mu_t * pow(M_inf, 0.4e1) * T_inf * T_inf * u_tau * u_tau * U * U * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(u_inf, -0.4e1) / Pr_t * pow(x, -0.2e1) / NU_SA / RHO / T / 0.196e3 + (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (Gamma - 0.1e1) * (u_eq_plus + y_plus * d_ueqplus_yplus) * cp * r_T * mu_t * kappa * M_inf * M_inf * T_inf * u_tau * u_tau * y * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / Pr_t * pow(x, -0.2e1) * pow(NU_SA, -0.2e1) / RHO / 0.196e3 - (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (-kappa * u_tau + 0.2e1 * alpha * y) * (Gamma - 0.1e1) * cp * r_T * mu_t * M_inf * M_inf * T_inf * y_plus * u_tau * d_ueqplus_yplus * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / Pr_t / y * pow(NU_SA, -0.2e1) / RHO + (-mu_t / Pr_t - mu / Pr) * (D2TDx2 + D2TDy2) * cp;

  // "Contribution  from the viscous/turbulence work to the total source term ----------------------------------------------"
  Scalar Q_E_work = -(0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (Gamma - 0.1e1) * r_T * mu_t * M_inf * M_inf * T_inf * y_plus * y_plus * u_tau * u_tau * d_ueqplus_yplus * d_ueqplus_yplus * U * U * (cos(A / u_inf * u_eq), Scalar(0.2e1)) * (u_inf, Scalar(-0.2e1)) * (y, Scalar(-0.2e1)) / NU_SA / RHO / T + (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (Gamma - 0.1e1) * (u_eq_plus + y_plus * d_ueqplus_yplus) * r_T * mu_t * M_inf * M_inf * T_inf * y_plus * u_tau * u_tau * d_ueqplus_yplus * U * V * (cos(A / u_inf * u_eq), Scalar(0.2e1)) * (u_inf, Scalar(-0.2e1)) / x / y / NU_SA / RHO / T / 0.42e2 - 0.4e1 / 0.3e1 * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (Gamma - 0.1e1) * r_T * mu_t * M_inf * M_inf * T_inf * y_plus * u_tau * d_ueqplus_yplus * U * V * V * cos(A / u_inf * u_eq) * (u_inf, Scalar(-0.2e1)) * (y, Scalar(-0.2e1)) / NU_SA / RHO / T - (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (Gamma - 0.1e1) * (u_eq_plus + y_plus * d_ueqplus_yplus, Scalar(0.2e1)) * r_T * mu_t * M_inf * M_inf * T_inf * u_tau * u_tau * U * U * (cos(A / u_inf * u_eq), Scalar(0.2e1)) * (x, Scalar(-0.2e1)) * (u_inf, Scalar(-0.2e1)) / NU_SA / RHO / T / 0.147e3 + (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (Gamma - 0.1e1) * (0.43e2 * y_plus * d_ueqplus_yplus - 0.2e1 * u_eq_plus) * r_T * mu_t * M_inf * M_inf * T_inf * u_tau * U * U * V * cos(A / u_inf * u_eq) * (u_inf, Scalar(-0.2e1)) / x / y / NU_SA / RHO / T / 0.42e2 - 0.15e2 / 0.196e3 * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (Gamma - 0.1e1) * (u_eq_plus + y_plus * d_ueqplus_yplus) * r_T * mu_t * M_inf * M_inf * T_inf * u_tau * U * V * V * cos(A / u_inf * u_eq) * (x, Scalar(-0.2e1)) * (u_inf, Scalar(-0.2e1)) / NU_SA / RHO / T - (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * mu_t * kappa * y_plus * u_tau * u_tau * d_ueqplus_yplus * U * cos(A / u_inf * u_eq) / y * (NU_SA, Scalar(-0.2e1)) / RHO + (0.8e1 * RHO * NU_SA - 0.6e1 * mu_t) * alpha * mu_t * y_plus * u_tau * d_ueqplus_yplus * U * cos(A / u_inf * u_eq) * (NU_SA, -0.2e1) / RHO + (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (y_plus * d_ueqplus_yplus - 0.2e1 * u_eq_plus) * mu_t * kappa * u_tau * u_tau * V * cos(A / u_inf * u_eq) / x * (NU_SA, -0.2e1) / RHO / 0.42e2 - (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (u_eq_plus + y_plus * d_ueqplus_yplus) * mu_t * kappa * u_tau * u_tau * y * U * cos(A / u_inf * u_eq) * (x, -0.2e1) * (NU_SA, -0.2e1) / RHO / 0.147e3 + 0.2e1 / 0.21e2 * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * alpha * (u_eq_plus + y_plus * d_ueqplus_yplus) * mu_t * u_tau * y * V * cos(A / u_inf * u_eq) / x * (NU_SA, -0.2e1) / RHO + 0.8e1 / 0.3e1 * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * alpha * mu_t * V * V * (NU_SA, -0.2e1) / RHO - (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (0.784e3 * x * x + 0.45e2 * y * y) * mu_t * kappa * u_tau * V * V * (x, -0.2e1) / y * (NU_SA, -0.2e1) / RHO / 0.588e3 - (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (-0.43e2 * kappa * u_tau + 0.90e2 * alpha * y) * mu_t * U * V / x * (NU_SA, -0.2e1) / RHO / 0.42e2 + (-0.2e1 * mu_t - 0.2e1 * mu) * (y_plus * y_plus * u_tau * u_tau * d_ueqplus_yplus * d_ueqplus_yplus * (cos(A / u_inf * u_eq), 0.2e1) * (y, -0.2e1) / 0.2e1 + (u_eq_plus + y_plus * d_ueqplus_yplus, 0.2e1) * u_tau * u_tau * (cos(A / u_inf * u_eq), 0.2e1) * (x, -0.2e1) / 0.294e3 - (0.43e2 * y_plus * d_ueqplus_yplus - 0.2e1 * u_eq_plus) * u_tau * V * cos(A / u_inf * u_eq) / x / y / 0.42e2 + (0.2e1 / 0.3e1 * D2uDx2 + D2vDxy / 0.6e1 + D2uDy2 / 0.2e1) * U + (D2uDxy / 0.6e1 + D2vDx2 / 0.2e1 + 0.2e1 / 0.3e1 * D2vDy2) * V + 0.225e3 / 0.392e3 * V * V * (x, Scalar(-0.2e1)) + 0.2e1 / 0.3e1 * V * V * (y, Scalar(-0.2e1)));

  // "Total source term ----------------------------------------------------------------------------------"
  Scalar Q_E = Q_E_work + Q_E_heat_flux + Q_E_convection;
  return(Q_E);

}



// ----------------------------------------
//   Manufactured Analytical Solutions
// ----------------------------------------

template <typename Scalar>
Scalar MASA::fans_sa_steady_wall_bounded<Scalar>::eval_exact_u(Scalar x,Scalar y)
{
  // this resets all constants set from parameters
  update(x,y);
  Scalar u_an;
  u_an = u_inf / A * sin(A / u_inf * u_eq);
  return u_an;
}

template <typename Scalar>
Scalar MASA::fans_sa_steady_wall_bounded<Scalar>::eval_exact_v(Scalar x,Scalar y)
{
  // this resets all constants set from parameters
  update(x,y);
  Scalar v_an;
  v_an = eta_v * u_tau * y / x / 0.14e2; 
  return v_an;
}

template <typename Scalar>
Scalar MASA::fans_sa_steady_wall_bounded<Scalar>::eval_exact_t(Scalar x,Scalar y)
{
  // this resets all constants set from parameters
  update(x,y);
  Scalar T_an;
  Scalar u_an = eval_exact_u(x,y);
  T_an = T_inf * (0.1e1 + r_T * (Gamma - 0.1e1) * M_inf * M_inf * (0.1e1 - u_an * u_an * pow(u_inf, Scalar(-0.2e1))) / 0.2e1);
  return T_an;
}

template <typename Scalar>
Scalar MASA::fans_sa_steady_wall_bounded<Scalar>::eval_exact_rho(Scalar x,Scalar y)
{
  // this resets all constants set from parameters
  update(x,y);
  Scalar rho_an;
  rho_an = p_0 / R / T;
  return rho_an;
}

template <typename Scalar>
Scalar MASA::fans_sa_steady_wall_bounded<Scalar>::eval_exact_nu(Scalar x,Scalar y)
{
  // this resets all constants set from parameters
  update(x,y);
  Scalar nu_an;
  nu_an = kappa * u_tau * y - alpha * y * y; 
  return nu_an;
}

template <typename Scalar>
Scalar MASA::fans_sa_steady_wall_bounded<Scalar>::eval_exact_p(Scalar x,Scalar y)
{
  // this resets all constants set from parameters
  update(x,y);
  return p_0;
}

// update function -- set all constants that are computed from parameters
template <typename Scalar>
Scalar MASA::fans_sa_steady_wall_bounded<Scalar>::update(Scalar x, Scalar y)
{

  // "---------------------------------------------"
  C1 = -0.1e1 / kappa * log(kappa) + C;
  u_inf = M_inf * sqrt(Gamma * R * T_inf);
  rho_inf = p_0 / R / T_inf;
  T_aw = T_inf * (0.1e1 + r_T * (Gamma - 0.1e1) * M_inf * M_inf / 0.2e1);
  rho_w = p_0 / R / T_aw;
  A = sqrt(0.1e1 - T_inf / T_aw);
  F_c = (T_aw / T_inf - 0.1e1) * pow(asin(A), -0.2e1);
  nu_w = mu / rho_w;

  // "---------------------------------------------"
  Re_x = rho_inf * u_inf * x / mu;
  c_f = C_cf / F_c * pow(0.1e1 / F_c * Re_x, -0.1e1 / 0.7e1);
  u_tau = u_inf * sqrt(c_f / 0.2e1);
  y_plus = y * u_tau / nu_w;
  u_eq_plus = 0.1e1 / kappa * log(0.1e1 + kappa * y_plus) + C1 * (0.1e1 - exp(-y_plus / eta1) - y_plus / eta1 * exp(-y_plus * b));
  u_eq = u_tau * u_eq_plus;
  U = u_inf / A * sin(A / u_inf * u_eq);
  V = eta_v * u_tau * y / x / 0.14e2;
  T = T_inf * (0.1e1 + r_T * (Gamma - 0.1e1) * M_inf * M_inf * (0.1e1 - U * U * pow(u_inf, -0.2e1)) / 0.2e1);
  RHO = p_0 / R / T;
  NU_SA = kappa * u_tau * y - alpha * y * y;
  chi = RHO * NU_SA / mu;
  f_v1 = pow(chi, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1));
  f_v2 = 0.1e1 - chi / (0.1e1 + chi * f_v1);
  mu_t = RHO * NU_SA * f_v1;

  // "---------------------------------------------"
  d_ueqplus_yplus = 0.1e1 / (0.1e1 + kappa * y_plus) + C1 * (exp(-y_plus / eta1) / eta1 - exp(-y_plus * b) / eta1 + y_plus * b * exp(-y_plus * b) / eta1);
  D2ueqDx2 = -u_tau * y_plus * y_plus * d_ueqplus_yplus / eta1 * pow(x, -0.2e1) / 0.196e3 + (0.17e2 * y_plus * d_ueqplus_yplus + 0.15e2 * u_eq_plus) * u_tau * pow(x, -0.2e1) / 0.196e3 - pow(x, -0.2e1) * ((b * b * eta1 * y_plus - 0.2e1 * b * eta1 - y_plus * b + 0.1e1) * C1 * u_tau * y_plus * y_plus * exp(-y_plus * b) * pow(eta1, -0.2e1) + (eta1 * kappa - kappa * y_plus - 0.1e1) * u_tau * y_plus * y_plus * pow(0.1e1 + kappa * y_plus, -0.2e1) / eta1) / 0.196e3;
  D2ueqDy2 = -u_tau * y_plus * y_plus * d_ueqplus_yplus / eta1 * pow(y, -0.2e1) - pow(y, -0.2e1) * ((b * b * eta1 * y_plus - 0.2e1 * b * eta1 - y_plus * b + 0.1e1) * C1 * u_tau * y_plus * y_plus * exp(-y_plus * b) * pow(eta1, -0.2e1) + (eta1 * kappa - kappa * y_plus - 0.1e1) * u_tau * y_plus * y_plus * pow(0.1e1 + kappa * y_plus, -0.2e1) / eta1);
  D2uDx2 = -u_tau * u_tau * A * A * pow(u_eq_plus + y_plus * d_ueqplus_yplus, 0.2e1) * U * pow(x, -0.2e1) * pow(u_inf, -0.2e1) / 0.196e3 + D2ueqDx2 * cos(A / u_inf * u_eq);
  D2uDy2 = -U * A * A * y_plus * y_plus * d_ueqplus_yplus * d_ueqplus_yplus * u_tau * u_tau * pow(u_inf, -0.2e1) * pow(y, -0.2e1) + cos(A / u_inf * u_eq) * D2ueqDy2;
  D2vDxy = -0.15e2 / 0.14e2 * V / x / y;
  D2TDx2 = -T_inf * r_T * (Gamma - 0.1e1) * M_inf * M_inf * pow(u_eq_plus + y_plus * d_ueqplus_yplus, 0.2e1) * u_tau * u_tau * (pow(cos(A / u_inf * u_eq), 0.2e1) - pow(sin(A / u_inf * u_eq), 0.2e1)) * pow(x, -0.2e1) * pow(u_inf, -0.2e1) / 0.196e3 - T_inf * r_T * (Gamma - 0.1e1) * M_inf * M_inf * D2ueqDx2 * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1);
  D2vDx2 = 0.435e3 / 0.196e3 * V * pow(x, -0.2e1);
  D2vDy2 = 0.0e0;
  D2TDy2 = -T_inf * r_T * (Gamma - 0.1e1) * M_inf * M_inf * y_plus * y_plus * d_ueqplus_yplus * d_ueqplus_yplus * u_tau * u_tau * (pow(cos(A / u_inf * u_eq), 0.2e1) - pow(sin(A / u_inf * u_eq), 0.2e1)) * pow(u_inf, -0.2e1) * pow(y, -0.2e1) - T_inf * r_T * (Gamma - 0.1e1) * M_inf * M_inf * D2ueqDy2 * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1);
  D2uDxy = A * A * (u_eq_plus + y_plus * d_ueqplus_yplus) * u_tau * u_tau * y_plus * d_ueqplus_yplus * U * pow(u_inf, -0.2e1) / x / y / 0.14e2 - u_tau * y_plus * d_ueqplus_yplus * cos(A / u_inf * u_eq) / x / y / 0.7e1 - y * cos(A / u_inf * u_eq) * D2ueqDy2 / x / 0.14e2;


  // "distance from wall"
  d = y;
  
  // magnitude of vorticity
  Omega = sqrt(pow(0.196e3 * x * x * y_plus * d_ueqplus_yplus * cos(A / u_inf * u_eq) + 0.15e2 * eta_v * y * y, 0.2e1) * u_tau * u_tau * pow(x, -0.4e1) * pow(y, -0.2e1)) / 0.196e3;


  // hacking in a few more constants
  c_w1 = c_b1 * pow(kappa, -0.2e1) + (0.1e1 + c_b2) / sigma;
  Sm_orig = NU_SA * pow(kappa, -0.2e1) * pow(d, -0.2e1) * f_v2;
  Sm1 = Sm_orig;
  Sm2 = Omega * (c_v2 * c_v2 * Omega + c_v3 * Sm_orig) / ((c_v3 + (-0.1e1) * 0.20e1 * c_v2) * Omega - Sm_orig);
  if (-c_v2 * Omega <= Sm_orig)
    Sm = Sm1;
  else
    Sm = Sm2;
  S_sa = Sm + Omega;
  cp = Gamma * R / (Gamma - 0.1e1);

  // dissipation stuff
  r = NU_SA / S_sa * pow(kappa, -0.2e1) * pow(d, -0.2e1);
  g = r + c_w2 * (pow(r, 0.6e1) - r);
  f_w = g * pow((0.1e1 + pow(c_w3, 0.6e1)) / (pow(g, 0.6e1) + pow(c_w3, 0.6e1)), 0.1e1 / 0.6e1);

  return 0;

}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::fans_sa_steady_wall_bounded);
MASA_INSTANTIATE_ALL(MASA::fans_sa_transient_free_shear);
