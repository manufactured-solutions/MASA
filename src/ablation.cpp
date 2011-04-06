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
// ablation.cpp: These are the MASA class member functions and constructors
//               For ablation+flow
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h> 

using namespace MASA;

/* ------------------------------------------------
 *   Ablation
 *   Navier Stokes 1D
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
MASA::navierstokes_ablation_1d_steady<Scalar>::navierstokes_ablation_1d_steady()
{
    this->mmsname = "navierstokes_ablation_1d_steady";
    this->dimension=1;

    this->register_var("R",&R);
    this->register_var("k",&k);

    this->register_var("u_0",&u_0);
    this->register_var("u_x",&u_x);
    this->register_var("a_ux",&a_ux);
    this->register_var("L",&L);
    this->register_var("Gamma",&Gamma);
    this->register_var("mu",&mu);

    this->register_var("T_0",&T_0);
    this->register_var("T_x",&T_x);
    this->register_var("a_Tx",&a_Tx);

    this->register_var("W_C",&W_C);
    this->register_var("W_C3",&W_C3);

    this->register_var("rho_N_0",&rho_N_0);
    this->register_var("rho_N_x",&rho_N_x);
    this->register_var("rho_N2_0",&rho_N2_0);
    this->register_var("rho_N2_x",&rho_N2_x);
    this->register_var("a_rho_N2_x",&a_rho_N2_x);

    this->register_var("rho_C3_0",&rho_C3_0);
    this->register_var("rho_C3_x",&rho_C3_x);
    this->register_var("a_rho_C3_x",&a_rho_C3_x);

    this->register_var("rho_C_0",&rho_C_0);
    this->register_var("rho_C_x",&rho_C_x);
    this->register_var("a_rho_C_x",&a_rho_C_x);

    this->register_var("rho_an_C3",&rho_an_C3);

    this->register_var("k_B",&k_B);
    this->register_var("beta_C3",&beta_C3);
    this->register_var("A_C3Enc",&A_C3Enc);
    this->register_var("D_C",&D_C);
    this->register_var("D_C3",&D_C3);
    this->register_var("m_C3",&m_C3);
    this->register_var("E_aC3nc",&E_aC3nc);

    this->register_var("qr",&qr);
    this->register_var("alpha",&alpha);
    this->register_var("a_rho_N_x",&a_rho_N_x);
    this->register_var("sigma",&sigma);
    this->register_var("epsilon",&epsilon);

    this->init_var();
  
}//done with constructor

template <typename Scalar>
int MASA::navierstokes_ablation_1d_steady<Scalar>::init_var()
{
  int err = 0;
  err += this->set_var("k",1.38);
  err += this->set_var("R",1.01);
  err += this->set_var("u_0",14.191);
  err += this->set_var("u_x",1.63);
  err += this->set_var("a_ux",.03);
  err += this->set_var("L",3.02);
  err += this->set_var("Gamma",16.1);
  err += this->set_var("mu",.091);

  err += this->set_var("T_0",12);
  err += this->set_var("T_x",12);
  err += this->set_var("a_Tx",12);

  err += this->set_var("W_C",12);
  err += this->set_var("W_C3",12);

  err += this->set_var("rho_N_0",12);
  err += this->set_var("rho_N_x",12);
  err += this->set_var("rho_N2_0",12);
  err += this->set_var("rho_N2_x",12);
  err += this->set_var("a_rho_N2_x",12);

  err += this->set_var("rho_C3_0",12);
  err += this->set_var("rho_C3_x",12);
  err += this->set_var("a_rho_C3_x",12);

  err += this->set_var("rho_C_0",12);
  err += this->set_var("rho_C_x",12);
  err += this->set_var("a_rho_C_x",12);

  err += this->set_var("rho_an_C3",12);

  err += this->set_var("k_B",12);
  err += this->set_var("beta_C3",12);
  err += this->set_var("A_C3Enc",12);
  err += this->set_var("D_C",12);
  err += this->set_var("D_C3",12);
  err += this->set_var("m_C3",12);
  err += this->set_var("E_aC3nc",12);

  err += this->set_var("qr",12);
  err += this->set_var("alpha",12);
  err += this->set_var("a_rho_N_x",12);
  err += this->set_var("sigma",12);
  err += this->set_var("epsilon",12);

  return err;

}

template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_q_rho_u(Scalar x)
{  
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
  Q_u = -a_rho_C3_x * PI * rho_C3_x * R * T * sin(a_rho_C3_x * PI * x / L) / L / W_C3 + a_rho_C_x * PI * rho_C_x * R * T * cos(a_rho_C_x * PI * x / L) / L / W_C + 0.4e1 / 0.3e1 * mu * a_ux * a_ux * PI * PI * u_x * sin(a_ux * PI * x / L) * pow(L, -0.2e1) + 0.2e1 * a_ux * PI * u_x * RHO * U * cos(a_ux * PI * x / L) / L;
  return(Q_u);

}

//flow energy equation
template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_q_rho_e(Scalar x)
{
  Scalar Q_e;
  Scalar RHO;
  Scalar RHO_N;
  Scalar RHO_N2;
  Scalar RHO_C;
  Scalar RHO_C3;
  Scalar U;
  Scalar T;
  Scalar P;
  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L);
  RHO_C = rho_C_0 + rho_C_x * sin(a_rho_C_x * PI * x / L);
  RHO_C3 = rho_C3_0 + rho_C3_x * cos(a_rho_C3_x * PI * x / L);
  RHO = RHO_N + RHO_N2;
  U = u_0 + u_x * sin(a_ux * PI * x / L);
  T = T_0 + T_x * cos(a_Tx * PI * x / L);
  P = R * T * (RHO_C / W_C + RHO_C3 / W_C3);
  Q_e = -0.4e1 / 0.3e1 * mu * a_ux * a_ux * PI * PI * u_x * u_x * pow(cos(a_ux * PI * x / L), 0.2e1) * pow(L, -0.2e1) - a_rho_C3_x * PI * rho_C3_x * R * T * U * sin(a_rho_C3_x * PI * x / L) / W_C3 / L + a_rho_C_x * PI * rho_C_x * R * T * U * cos(a_rho_C_x * PI * x / L) / W_C / L - a_Tx * PI * T_x * R * RHO_C3 * U * sin(a_Tx * PI * x / L) / W_C3 / L - a_Tx * PI * T_x * R * RHO_C * U * sin(a_Tx * PI * x / L) / W_C / L - a_Tx * PI * T_x * R * RHO * U * sin(a_Tx * PI * x / L) / (Gamma - 0.1e1) / L + 0.4e1 / 0.3e1 * mu * a_ux * a_ux * PI * PI * u_x * U * sin(a_ux * PI * x / L) * pow(L, -0.2e1) + a_ux * PI * u_x * R * RHO * T * cos(a_ux * PI * x / L) / (Gamma - 0.1e1) / L + 0.3e1 / 0.2e1 * a_ux * PI * u_x * RHO * U * U * cos(a_ux * PI * x / L) / L + k * a_Tx * a_Tx * PI * PI * T_x * cos(a_Tx * PI * x / L) * pow(L, -0.2e1) + a_ux * PI * u_x * P * cos(a_ux * PI * x / L) / L - (a_rho_C3_x * rho_C3_x * sin(a_rho_C3_x * PI * x / L) - a_rho_C_x * rho_C_x * cos(a_rho_C_x * PI * x / L)) * PI * R * T * U / (Gamma - 0.1e1) / L - (a_rho_C3_x * rho_C3_x * sin(a_rho_C3_x * PI * x / L) - a_rho_C_x * rho_C_x * cos(a_rho_C_x * PI * x / L)) * PI * pow(U, 0.3e1) / L / 0.2e1;
  return(Q_e);
}

// flow equation for Carbon
template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_q_rho_C(Scalar x)
{
  Scalar Q_rho_C;
  Scalar U;
  U = u_0 + u_x * sin(a_ux * PI * x / L);
  Q_rho_C = a_rho_C_x * PI * rho_C_x * U * cos(a_rho_C_x * PI * x / L) / L;
  return(Q_rho_C);
}

// flow equation for Carbon3
template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_q_rho_C3(Scalar x)
{
  Scalar Q_rho_C3;
  Scalar U;
  U = u_0 + u_x * sin(a_ux * PI * x / L);
  Q_rho_C3 = -a_rho_C3_x * PI * rho_C3_x * U * sin(a_rho_C3_x * PI * x / L) / L;
  return(Q_rho_C3);
}

// this is the ablation energy equation
template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_q_e(Scalar x,Scalar (*in_func)(Scalar))
{
  Scalar Q_e;
  Scalar P;
  Scalar T;
  Scalar RHO;
  Scalar RHO_C;
  Scalar RHO_C3;
  Scalar MF_C3;
  Scalar MF_C3E;
  Scalar Mdot_C3C;
  Scalar H_C3;
  RHO_C3 = rho_C3_0 + rho_C3_x * cos(a_rho_C3_x * PI * x / L);
  RHO_C = rho_C_0 + rho_C_x * sin(a_rho_C_x * PI * x / L);
  RHO = RHO_C + RHO_C3;
  T = T_0 + T_x * cos(a_Tx * PI * x / L);
  P = R * T * (RHO_C / W_C + RHO_C3 / W_C3);
  MF_C3 = RHO_C3 / RHO;
  MF_C3E = A_C3Enc * exp(-E_aC3nc / T) / P;
  Mdot_C3C = sqrt(T * k_B / PI / m_C3) * sqrt(0.2e1) * (-MF_C3 + MF_C3E) * RHO * beta_C3 / 0.2e1;
  //H_C3 = Function_to_Calculate_h_C3;
  H_C3 = in_func(T);

  Q_e = k * a_Tx * PI * T_x * sin(a_Tx * PI * x / L) / L - sigma * epsilon * pow(T, 0.4e1) - Mdot_C3C * H_C3 + alpha * qr;
  return(Q_e);
}

// ablation equation for Carbon
template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_q_C(Scalar x)
{
  int Q_rho_C;
  Q_rho_C = 0*x;
  return(Q_rho_C);
}

// ablation equation for Carbon3
template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_q_C3(Scalar x)
{
  Scalar Q_rho_C3;
  Scalar RHO;
  Scalar RHO_C;
  Scalar RHO_C3;
  Scalar T;
  Scalar P;
  Scalar MF_C3;
  Scalar MF_C3E;
  Scalar Mdot_C3C;
  RHO_C3 = rho_C3_0 + rho_C3_x * cos(a_rho_C3_x * PI * x / L);
  RHO_C = rho_C_0 + rho_C_x * sin(a_rho_C_x * PI * x / L);
  RHO = RHO_C + RHO_C3;
  T = T_0 + T_x * cos(a_Tx * PI * x / L);
  P = R * T * (RHO_C / W_C + RHO_C3 / W_C3);
  MF_C3 = RHO_C3 / RHO;
  MF_C3E = A_C3Enc * exp(-E_aC3nc / T) / P;
  Mdot_C3C = sqrt(T * k_B / PI / m_C3) * sqrt(0.2e1) * (-MF_C3 + MF_C3E) * RHO * beta_C3 / 0.2e1;
  Q_rho_C3 = a_rho_C_x * PI * rho_C_x * D_C3 * RHO_C3 * cos(a_rho_C_x * PI * x / L) / L / RHO + a_rho_C3_x * PI * rho_C3_x * D_C3 * RHO_C * sin(a_rho_C3_x * PI * x / L) / L / RHO + a_rho_C_x * (D_C - D_C3) * PI * rho_C_x * RHO_C3 * RHO_C3 * cos(a_rho_C_x * PI * x / L) / L * pow(RHO, -0.2e1) + a_rho_C3_x * (D_C - D_C3) * PI * rho_C3_x * RHO_C * RHO_C3 * sin(a_rho_C3_x * PI * x / L) / L * pow(RHO, -0.2e1) + Mdot_C3C * (MF_C3 - 0.1e1);
  return(Q_rho_C3);
}

template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_q_u_boundary(Scalar x)
{
  Scalar Q_u_boundary;
  Scalar RHO;
  Scalar RHO_C;
  Scalar RHO_C3;
  Scalar U;
  Scalar T;
  Scalar P;
  Scalar MF_C3;
  Scalar MF_C3E;
  Scalar Mdot_C3C;
  RHO_C3 = rho_C3_0 + rho_C3_x * cos(a_rho_C3_x * PI * x / L);
  RHO_C = rho_C_0 + rho_C_x * sin(a_rho_C_x * PI * x / L);
  RHO = RHO_C + RHO_C3;
  //nick adding U here:
  U = u_0 + u_x * sin(a_ux * PI * x / L);
  T = T_0 + T_x * cos(a_Tx * PI * x / L);
  P = R * T * (RHO_C / W_C + RHO_C3 / W_C3);
  MF_C3 = RHO_C3 / RHO;
  MF_C3E = A_C3Enc * exp(-E_aC3nc / T) / P;
  Mdot_C3C = sqrt(T * k_B / PI / m_C3) * sqrt(0.2e1) * (-MF_C3 + MF_C3E) * RHO * beta_C3 / 0.2e1;
  Q_u_boundary = -Mdot_C3C / RHO + U;
  return(Q_u_boundary);
}

/* ------------------------------------------------
 * 
 *
 *    Manufactured Analytical Solutions
 * 
 * -----------------------------------------------
 */ 

template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_exact_u(Scalar x)
{
  Scalar u_an = u_0 + u_x * sin(a_ux * pi * x / L);
  return u_an;
}

template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_exact_t(Scalar x)
{
  Scalar T_an = T_0 + T_x * cos(a_Tx * pi * x / L);
  return T_an;
}

template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_exact_rho(Scalar x)
{
  Scalar rho_an = rho_C_0 + rho_C_x * sin(a_rho_C_x * pi * x / L) + rho_C3_0 + rho_C3_x * cos(a_rho_C3_x * pi * x / L);
  return rho_an;
}

template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_exact_rho_C(Scalar x)
{
  Scalar rho_an_C = rho_C_0 + rho_C_x * sin(a_rho_C_x * pi * x / L);
  return rho_an_C;
}

template <typename Scalar>
Scalar MASA::navierstokes_ablation_1d_steady<Scalar>::eval_exact_rho_C3(Scalar x)
{
  Scalar rho_an_C3 = rho_C3_0 + rho_C3_x * cos(a_rho_C3_x * pi * x / L);
  return rho_an_C3;
}

MASA_INSTANTIATE_ALL(MASA::navierstokes_ablation_1d_steady);
