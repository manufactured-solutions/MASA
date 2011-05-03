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
// $Author: 
// $Id:
//
// euler_chem.cpp: 
//          These are the MASA class member functions and constructors
//          For the Euler Equations + Chemistry
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h> 

using namespace MASA;
/* ------------------------------------------------
 *
 *         EULER EQUATION 1D
 *
 *
 *
 * -----------------------------------------------
 */ 
template <typename Scalar>
MASA::euler_chem_1d<Scalar>::euler_chem_1d()
{
  this->mmsname = "euler_chem_1d";
  this->dimension=1;

  this->register_var("R",&R);
  this->register_var("Cf1_N",&Cf1_N);
  this->register_var("Cf1_N2",&Cf1_N2);

  this->register_var("etaf1_N",&etaf1_N);
  this->register_var("etaf1_N2",&etaf1_N2);

  this->register_var("Ea_N",&Ea_N);
  this->register_var("Ea_N2",&Ea_N2);

  this->register_var("R_N",&R_N);
  this->register_var("R_N2",&R_N2);

  this->register_var("theta_v_N2",&theta_v_N2);
  this->register_var("M_N",&M_N);

  this->register_var("h0_N",&h0_N);
  this->register_var("h0_N2",&h0_N2);
  this->register_var("L",&L);

  this->register_var("u_0",&u_0);
  this->register_var("u_x",&u_x);
  this->register_var("a_ux",&a_ux);

  this->register_var("rho_N_0",&rho_N_0);
  this->register_var("rho_N_x",&rho_N_x);
  this->register_var("a_rho_N_x",&a_rho_N_x);

  this->register_var("rho_N2_0",&rho_N2_0);
  this->register_var("rho_N2_x",&rho_N2_x);
  this->register_var("a_rho_N2_x",&a_rho_N2_x);

  this->register_var("T_0",&T_0);
  this->register_var("T_x",&T_x);
  this->register_var("a_Tx",&a_Tx);

  this->init_var();
}

template <typename Scalar>
int MASA::euler_chem_1d<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("R",287);
  err += this->set_var("Cf1_N",1.01);
  err += this->set_var("Cf1_N2",3.11);

  err += this->set_var("etaf1_N",0.21);
  err += this->set_var("etaf1_N2",1.01);

  err += this->set_var("Ea_N",0.21);
  err += this->set_var("Ea_N2",2.71);

  err += this->set_var("R_N",0.60);
  err += this->set_var("R_N2",0.40);

  err += this->set_var("theta_v_N2",1.01);
  err += this->set_var("M_N",5.01);

  err += this->set_var("h0_N",1.31);
  err += this->set_var("h0_N2",2.21);
  err += this->set_var("L",10.00);

  err += this->set_var("u_0",70);
  err += this->set_var("u_x", 4);
  err += this->set_var("a_ux",5/3);

  err += this->set_var("rho_N_0",0.81);
  err += this->set_var("rho_N_x",0.11);
  err += this->set_var("a_rho_N_x",8.01);

  err += this->set_var("rho_N2_0",4.11);
  err += this->set_var("rho_N2_x",0.41);
  err += this->set_var("a_rho_N2_x",10.01);

  err += this->set_var("T_0",10.01);
  err += this->set_var("T_x",7/4);
  err += this->set_var("a_Tx",1.21);

  return err;
}


// ----------------------------------------
//   Source Terms
// ----------------------------------------

template <typename Scalar>
Scalar MASA::euler_chem_1d<Scalar>::eval_q_rho_u(Scalar x)
{
  //Scalar x, Scalar R_N

  Scalar Q_u;
  Scalar RHO;
  Scalar RHO_N;
  Scalar RHO_N2;
  Scalar U;
  Scalar T;

  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L);
  RHO = RHO_N + RHO_N2;
  U = u_0 + u_x * sin(a_ux * PI * x / L);
  T = T_0 + T_x * cos(a_Tx * PI * x / L);

  Q_u = a_rho_N_x * PI * rho_N_x * U * U * cos(a_rho_N_x * PI * x / L) / L - a_rho_N2_x * PI * rho_N2_x * U * U * sin(a_rho_N2_x * PI * x / L) / L - a_Tx * PI * T_x * R_N * RHO_N * sin(a_Tx * PI * x / L) / L - a_Tx * PI * T_x * R_N * RHO_N2 * sin(a_Tx * PI * x / L) / L / Scalar(0.2e1) + Scalar(0.2e1) * a_ux * PI * u_x * RHO * U * cos(a_ux * PI * x / L) / L - (Scalar(-0.2e1) * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) + a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R_N * T / L / Scalar(0.2e1);
  return(Q_u);
}

template <typename Scalar>
Scalar MASA::euler_chem_1d<Scalar>::eval_q_rho_e(Scalar x)
{
  /*
  double x,
  double R_N,
  double R_N2,
  double h0_N,
  double h0_N2,
  double theta_v_N2)
  */

  Scalar Q_e;
  Scalar RHO;
  Scalar RHO_N;
  Scalar RHO_N2;
  Scalar U;
  Scalar T;
  Scalar alpha;
  Scalar E_vib_N2;
  Scalar C  = -0.2e1;
  Scalar C2 =  0.3e1;

  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L);
  RHO = RHO_N + RHO_N2;
  U = u_0 + u_x * sin(a_ux * PI * x / L);
  T = T_0 + T_x * cos(a_Tx * PI * x / L);
  alpha = exp(theta_v_N2 / T);
  E_vib_N2 = R_N2 * theta_v_N2 / (alpha - 0.1e1);

  Q_e = -a_Tx * PI * T_x * alpha * theta_v_N2 * E_vib_N2 * RHO_N2 * U * sin(a_Tx * PI * x / L) / L / (alpha - Scalar(0.1e1)) * pow(T, C) - Scalar(0.5e1) / Scalar(0.2e1) * a_Tx * PI * T_x * R_N * RHO_N * U * sin(a_Tx * PI * x / L) / L - Scalar(0.7e1) / Scalar(0.4e1) * a_Tx * PI * T_x * R_N * RHO_N2 * U * sin(a_Tx * PI * x / L) / L + Scalar(0.5e1) / Scalar(0.2e1) * a_ux * PI * u_x * R_N * RHO_N * T * cos(a_ux * PI * x / L) / L + Scalar(0.7e1) / Scalar(0.4e1) * a_ux * PI * u_x * R_N * RHO_N2 * T * cos(a_ux * PI * x / L) / L + Scalar(0.3e1) / Scalar(0.2e1) * a_ux * PI * u_x * RHO * U * U * cos(a_ux * PI * x / L) / L - a_rho_N2_x * PI * rho_N2_x * E_vib_N2 * U * sin(a_rho_N2_x * PI * x / L) / L + a_ux * PI * u_x * E_vib_N2 * RHO_N2 * cos(a_ux * PI * x / L) / L + (h0_N * RHO_N + h0_N2 * RHO_N2) * a_ux * PI * u_x * cos(a_ux * PI * x / L) / L - (Scalar(-0.10e2) * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) + Scalar(0.7e1) * a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R_N * U * T / L / Scalar(0.4e1) - (-a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) + a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * pow(U, C2) / L / Scalar(0.2e1) - (-a_rho_N_x * rho_N_x * h0_N * cos(a_rho_N_x * PI * x / L) + a_rho_N2_x * rho_N2_x * h0_N2 * sin(a_rho_N2_x * PI * x / L)) * PI * U / L;

  return(Q_e);
}

template <typename Scalar>
Scalar MASA::euler_chem_1d<Scalar>::eval_q_rho_N(Scalar x, Scalar (*in_func)(Scalar))
{

  Scalar Q_rho_N;
  Scalar RHO_N;
  Scalar RHO_N2;
  Scalar U;
  Scalar T;
  Scalar kf1_N;
  Scalar kf1_N2;
  Scalar K_eq;
  Scalar C  = -0.2e1;

  // Calculate Equilibrium Konstant
  T = T_0 + T_x * cos(a_Tx * PI * x / L);
  K_eq = in_func(T);

  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L);
  kf1_N = Cf1_N * pow(T, etaf1_N) * exp(-Ea_N / R / T);
  kf1_N2 = Cf1_N2 * pow(T, etaf1_N2) * exp(-Ea_N2 / R / T);

  Q_rho_N = a_rho_N_x * PI * rho_N_x * U * cos(a_rho_N_x * PI * x / L) / L + a_ux * PI * u_x * RHO_N * cos(a_ux * PI * x / L) / L + (Scalar(0.2e1) * kf1_N * RHO_N + kf1_N2 * RHO_N2) * RHO_N * RHO_N * pow(M_N, C) / K_eq - (Scalar(0.2e1) * kf1_N * RHO_N + kf1_N2 * RHO_N2) * RHO_N2 / M_N / Scalar(0.2e1);

  return(Q_rho_N);
}

template <typename Scalar>
Scalar MASA::euler_chem_1d<Scalar>::eval_q_rho_N2(Scalar x, Scalar (*in_func)(Scalar))
{

  Scalar Q_rho_N2;
  Scalar RHO_N;
  Scalar RHO_N2;
  Scalar U;
  Scalar T;
  Scalar kf1_N;
  Scalar kf1_N2;
  Scalar K_eq;
  Scalar C  = -0.2e1;

  // Calculate Equilibrium Konstant
  T = T_0 + T_x * cos(a_Tx * PI * x / L);
  K_eq = in_func(T);
  
  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L);
  kf1_N = Cf1_N * pow(T, etaf1_N) * exp(-Ea_N / R / T);
  kf1_N2 = Cf1_N2 * pow(T, etaf1_N2) * exp(-Ea_N2 / R / T);

  Q_rho_N2 = -a_rho_N2_x * PI * rho_N2_x * U * sin(a_rho_N2_x * PI * x / L) / L + a_ux * PI * u_x * RHO_N2 * cos(a_ux * PI * x / L) / L - (Scalar(0.2e1) * kf1_N * RHO_N + kf1_N2 * RHO_N2) * RHO_N * RHO_N * pow(M_N, C) / K_eq + (Scalar(0.2e1) * kf1_N * RHO_N + kf1_N2 * RHO_N2) * RHO_N2 / M_N / Scalar(0.2e1);
  return(Q_rho_N2);
}

// ----------------------------------------
//   Analytical Terms
// ----------------------------------------

template <typename Scalar>
Scalar MASA::euler_chem_1d<Scalar>::eval_exact_t(Scalar x)
{
  Scalar T_an = T_0 + T_x * cos(a_Tx * pi * x / L);
  return T_an;
}

template <typename Scalar>
Scalar MASA::euler_chem_1d<Scalar>::eval_exact_u(Scalar x)
{
  Scalar u_an = u_0 + u_x * sin(a_ux * pi * x / L);
  return u_an;
}

template <typename Scalar>
Scalar MASA::euler_chem_1d<Scalar>::eval_exact_rho(Scalar x)
{
  Scalar rho_an = rho_N_0 + rho_N_x * sin(a_rho_N_x * pi * x / L) + rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * pi * x / L);
  return rho_an;
}

template <typename Scalar>
Scalar MASA::euler_chem_1d<Scalar>::eval_exact_rho_N(Scalar x)
{
  Scalar rho_an_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * pi * x / L);
  return rho_an_N;
}

template <typename Scalar>
Scalar MASA::euler_chem_1d<Scalar>::eval_exact_rho_N2(Scalar x)
{
  Scalar rho_an_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * pi * x / L);
  return rho_an_N2;
}

// ----------------------------------------
//   Gradients
// ----------------------------------------

/*
  grad_rho_an_N[0] = rho_N_x * cos(a_rho_N_x * pi * x / L) * a_rho_N_x * pi / L;
  grad_rho_an_N[1] = 0;
  grad_rho_an_N[2] = 0;
  grad_rho_an_N2[0] = -rho_N2_x * sin(a_rho_N2_x * pi * x / L) * a_rho_N2_x * pi / L;
  grad_rho_an_N2[1] = 0;
  grad_rho_an_N2[2] = 0;
  grad_rho_an[0] = rho_N_x * cos(a_rho_N_x * pi * x / L) * a_rho_N_x * pi / L - rho_N2_x * sin(a_rho_N2_x * pi * x / L) * a_rho_N2_x * pi / L;
  grad_rho_an[1] = 0;
  grad_rho_an[2] = 0;
  grad_u_an[0] = u_x * cos(a_ux * pi * x / L) * a_ux * pi / L;
  grad_u_an[1] = 0;
  grad_u_an[2] = 0;
  grad_T_an[0] = -T_x * sin(a_Tx * pi * x / L) * a_Tx * pi / L;
  grad_T_an[1] = 0;
  grad_T_an[2] = 0;
*/

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::euler_chem_1d);

