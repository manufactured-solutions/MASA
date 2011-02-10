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
// $Author:
// $Id:
//
// euler_chem.cpp: program that tests euler with chemistry
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <math.h>
#include <limits>
#include <iostream>
#include <stdlib.h>

using namespace MASA;
using namespace std;

template<typename Scalar>
Scalar nancheck(Scalar x)
{
  if(isnan(x))
    {
      cout << "MASA REGRESSION FAILURE:: nan found!\n";
      exit(1);
    }
  return 1;
}


template <typename Scalar>
Scalar SourceQ_rho_u(Scalar x, Scalar R_N)
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

  Q_u = a_rho_N_x * PI * rho_N_x * U * U * cos(a_rho_N_x * PI * x / L) / L - a_rho_N2_x * PI * rho_N2_x * U * U * sin(a_rho_N2_x * PI * x / L) / L - a_Tx * PI * T_x * R_N * RHO_N * sin(a_Tx * PI * x / L) / L - a_Tx * PI * T_x * R_N * RHO_N2 * sin(a_Tx * PI * x / L) / L / 0.2e1 + 0.2e1 * a_ux * PI * u_x * RHO * U * cos(a_ux * PI * x / L) / L - (-0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) + a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R_N * T / L / 0.2e1;
  return(Q_u);
}

template <typename Scalar>
Scalar SourceQ_rho_e(Scalar x,
		    Scalar R_N,
		    Scalar R_N2,
		    Scalar h0_N,
		    Scalar h0_N2,
		    Scalar theta_v_N2)
{
  Scalar Q_e;
  Scalar RHO;
  Scalar RHO_N;
  Scalar RHO_N2;
  Scalar U;
  Scalar T;
  Scalar alpha;
  Scalar E_vib_N2;

  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L);
  RHO = RHO_N + RHO_N2;
  U = u_0 + u_x * sin(a_ux * PI * x / L);
  T = T_0 + T_x * cos(a_Tx * PI * x / L);
  alpha = exp(theta_v_N2 / T);
  E_vib_N2 = R_N2 * theta_v_N2 / (alpha - 0.1e1);

  Q_e = -a_Tx * PI * T_x * alpha * theta_v_N2 * E_vib_N2 * RHO_N2 * U * sin(a_Tx * PI * x / L) / L / (alpha - 0.1e1) * pow(T, -0.2e1) - 0.5e1 / 0.2e1 * a_Tx * PI * T_x * R_N * RHO_N * U * sin(a_Tx * PI * x / L) / L - 0.7e1 / 0.4e1 * a_Tx * PI * T_x * R_N * RHO_N2 * U * sin(a_Tx * PI * x / L) / L + 0.5e1 / 0.2e1 * a_ux * PI * u_x * R_N * RHO_N * T * cos(a_ux * PI * x / L) / L + 0.7e1 / 0.4e1 * a_ux * PI * u_x * R_N * RHO_N2 * T * cos(a_ux * PI * x / L) / L + 0.3e1 / 0.2e1 * a_ux * PI * u_x * RHO * U * U * cos(a_ux * PI * x / L) / L - a_rho_N2_x * PI * rho_N2_x * E_vib_N2 * U * sin(a_rho_N2_x * PI * x / L) / L + a_ux * PI * u_x * E_vib_N2 * RHO_N2 * cos(a_ux * PI * x / L) / L + (h0_N * RHO_N + h0_N2 * RHO_N2) * a_ux * PI * u_x * cos(a_ux * PI * x / L) / L - (-0.10e2 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) + 0.7e1 * a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R_N * U * T / L / 0.4e1 - (-a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) + a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * pow(U, 0.3e1) / L / 0.2e1 - (-a_rho_N_x * rho_N_x * h0_N * cos(a_rho_N_x * PI * x / L) + a_rho_N2_x * rho_N2_x * h0_N2 * sin(a_rho_N2_x * PI * x / L)) * PI * U / L;

  return(Q_e);
}

template <typename Scalar>
Scalar SourceQ_rho_N(Scalar x,
		    Scalar M_N,
		    Scalar h0_N,
		    Scalar h0_N2,
		    Scalar Cf1_N,
		    Scalar Cf1_N2,
		    Scalar etaf1_N,
		    Scalar etaf1_N2,
		    Scalar Ea_N,
		    Scalar Ea_N2,
		    Scalar Function_to_Calculate_K)
{
  Scalar Q_rho_N;
  Scalar RHO_N;
  Scalar RHO_N2;
  Scalar U;
  Scalar T;
  Scalar kf1_N;
  Scalar kf1_N2;
  Scalar K_eq;

  K_eq = Function_to_Calculate_K;
  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L);
  T = T_0 + T_x * cos(a_Tx * PI * x / L);
  kf1_N = Cf1_N * pow(T, etaf1_N) * exp(-Ea_N / R / T);
  kf1_N2 = Cf1_N2 * pow(T, etaf1_N2) * exp(-Ea_N2 / R / T);

  Q_rho_N = a_rho_N_x * PI * rho_N_x * U * cos(a_rho_N_x * PI * x / L) / L + a_ux * PI * u_x * RHO_N * cos(a_ux * PI * x / L) / L + (0.2e1 * kf1_N * RHO_N + kf1_N2 * RHO_N2) * RHO_N * RHO_N * pow(M_N, -0.2e1) / K_eq - (0.2e1 * kf1_N * RHO_N + kf1_N2 * RHO_N2) * RHO_N2 / M_N / 0.2e1;

  return(Q_rho_N);
}

template <typename Scalar>
Scalar SourceQ_rho_N2(Scalar x,
		     Scalar M_N,
		     Scalar h0_N,
		     Scalar h0_N2,
		     Scalar K,
		     Scalar Cf1_N,
		     Scalar Cf1_N2,
		     Scalar etaf1_N,
		     Scalar etaf1_N2,
		     Scalar Ea_N,
		     Scalar Ea_N2,
		     Scalar Function_to_Calculate_K)
{
  
  Scalar Q_rho_N2;
  Scalar RHO_N;
  Scalar RHO_N2;
  Scalar U;
  Scalar T;
  Scalar kf1_N;
  Scalar kf1_N2;
  Scalar K_eq;

  K_eq = Function_to_Calculate_K;
  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L);
  T = T_0 + T_x * cos(a_Tx * PI * x / L);
  kf1_N = Cf1_N * pow(T, etaf1_N) * exp(-Ea_N / R / T);
  kf1_N2 = Cf1_N2 * pow(T, etaf1_N2) * exp(-Ea_N2 / R / T);

  Q_rho_N2 = -a_rho_N2_x * PI * rho_N2_x * U * sin(a_rho_N2_x * PI * x / L) / L + a_ux * PI * u_x * RHO_N2 * cos(a_ux * PI * x / L) / L - (0.2e1 * kf1_N * RHO_N + kf1_N2 * RHO_N2) * RHO_N * RHO_N * pow(M_N, -0.2e1) / K_eq + (0.2e1 * kf1_N * RHO_N + kf1_N2 * RHO_N2) * RHO_N2 / M_N / 0.2e1;
  return(Q_rho_N2);
}

// ----------------------------------------
//   Analytical Terms
// ----------------------------------------

template <typename Scalar>
Scalar anQ_t(Scalar x)
{
  Scalar T_an = T_0 + T_x * cos(a_Tx * pi * x / L);
  return T_an;
}

template <typename Scalar>
Scalar anQ_u(Scalar x)
{
  Scalar u_an = u_0 + u_x * sin(a_ux * pi * x / L);
  return u_an;
}

template <typename Scalar>
Scalar anQ_rho(Scalar x)
{
  Scalar rho_an = rho_N_0 + rho_N_x * sin(a_rho_N_x * pi * x / L) + rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * pi * x / L);
  return rho_an;
}

template <typename Scalar>
Scalar anQ_rho_N(Scalar x)
{
  Scalar rho_an_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * pi * x / L);
  return rho_an_N;
}

template <typename Scalar>
Scalar anQ_rho_N2(Scalar x)
{
  Scalar rho_an_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * pi * x / L);
  return rho_an_N2;
}

// ----------------------------------------
//   Regresssion
// ----------------------------------------

template<typename Scalar>
Scalar temp_function(Scalar T)
{
  // hackish functional here
  // That's an eyeballed fit (focusing on the 5000K-6000K range) 
  // for the equilibrium constant for N2->N+N dissociation
  Scalar K = exp(4+(T-6000)/500);
  return K;
}

template<typename Scalar>
int run_regression()
{

  Scalar threshold = 5 * numeric_limits<Scalar>::epsilon();

  //variables 
  Scalar u_0;
  Scalar u_x;
  Scalar rho_0;
  Scalar rho_x;
  Scalar p_0;
  Scalar p_x;
  Scalar a_px;
  Scalar a_rhox;
  Scalar a_ux;
  Scalar Gamma;
  Scalar mu;
  Scalar L;

  // parameters
  Scalar x;

  //problem size
  int nx = 200;  // number of points
  int lx=10;     // length
  Scalar dx=Scalar(lx)/Scalar(nx);

  // solutions
  Scalar ufield,ufield2,ufield3;
  Scalar efield,efield2,efield3;
  Scalar rho,rho2,rho3;
  Scalar gradx,grady,gradz,gradp,gradrho;

  Scalar exact_u,exact_u2,exact_u3;
  Scalar exact_p,exact_p2,exact_p3;
  Scalar exact_rho,exact_rho2,exact_rho3;

  // initalize
  masa_init<Scalar>("euler-test","euler_chem_1d");

  // initialize the default parameters
  masa_init_param<Scalar>();

  // check that all terms have been initialized
  masa_sanity_check<Scalar>();

  return 0;

} // done with tests



int main()
{
  int err=0;

  err += run_regression<double>();
  //err += run_regression<long double>();

  return err;
}
