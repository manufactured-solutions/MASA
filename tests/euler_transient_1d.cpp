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

// source terms

template <typename Scalar>
Scalar eval_q_rho_u(Scalar x,Scalar t)
{
  double Q_u_t;
  double RHO;
  double U;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_t * sin(a_rhot * pi * t / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_t * cos(a_ut * pi * t / L);

  Q_u_t = cos(a_rhox * PI * x / L) * a_rhox * PI * rho_x * U * U / L + 0.2e1 * cos(a_ux * PI * x / L) * RHO * a_ux * PI * u_x * U / L + cos(a_rhot * PI * t / L) * a_rhot * PI * rho_t * U / L - sin(a_ut * PI * t / L) * a_ut * PI * u_t * RHO / L - sin(a_px * PI * x / L) * a_px * PI * p_x / L;
  return(Q_u_t);


}


template <typename Scalar>
Scalar eval_q_rho_e(Scalar x,Scalar t)
{
  double Q_e_t;
  double RHO;
  double U;
  double P;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_t * sin(a_rhot * pi * t / L);
  P = p_0 + p_x * cos(a_px * pi * x / L) + p_t * cos(a_pt * pi * t / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_t * cos(a_ut * pi * t / L);

  Q_e_t = cos(a_rhot * PI * t / L) * a_rhot * PI * rho_t * U * U / L / 0.2e1 - sin(a_ut * PI * t / L) * a_ut * PI * u_t * RHO * U / L - sin(a_pt * PI * t / L) * a_pt * PI * p_t / (Gamma - 0.1e1) / L + cos(a_rhox * PI * x / L) * pow(U, 0.3e1) * a_rhox * PI * rho_x / L / 0.2e1 + cos(a_ux * PI * x / L) * P * a_ux * PI * u_x * Gamma / (Gamma - 0.1e1) / L + 0.3e1 / 0.2e1 * cos(a_ux * PI * x / L) * RHO * U * U * a_ux * PI * u_x / L - sin(a_px * PI * x / L) * U * a_px * PI * p_x * Gamma / (Gamma - 0.1e1) / L;
  return(Q_e_t);

}

template <typename Scalar>
Scalar eval_q_rho(Scalar x,Scalar t)
{

  double Q_rho_t;
  double RHO;
  double U;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_t * sin(a_rhot * pi * t / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_t * cos(a_ut * pi * t / L);

  Q_rho_t = cos(a_rhot * PI * t / L) * a_rhot * PI * rho_t / L + cos(a_ux * PI * x / L) * RHO * a_ux * PI * u_x / L + cos(a_rhox * PI * x / L) * U * a_rhox * PI * rho_x / L;
  return(Q_rho_t);

}

/* ------------------------------------------------
 *
 *
 *   Analytical terms
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
Scalar eval_exact_u(Scalar x,Scalar t)
{
  Scalar exact_u;
  exact_u = u_0 + u_x * sin(a_ux * pi * x / L) + u_t * cos(a_ut * pi * t / L);
  return exact_u;
}

template <typename Scalar>
Scalar eval_exact_p(Scalar x,Scalar t)
{
  Scalar exact_p;
  exact_p = p_0 + p_x * cos(a_px * pi * x / L) + p_t * cos(a_pt * pi * t / L);
  return exact_p;
}

template <typename Scalar>
Scalar eval_exact_rho(Scalar x,Scalar t)
{
  Scalar exact_rho;
  exact_rho = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_t * sin(a_rhot * pi * t / L);
  return exact_rho;
}

template<typename Scalar>
int run_regression()
{

  Scalar threshold = 5 * numeric_limits<Scalar>::epsilon();

  Scalar ufield,ufield2,ufield3;
  Scalar efield,efield2,efield3;
  Scalar rhofield,rhofield2,rhofield3;
  
  Scalar exact_u;
  Scalar exact_p;
  Scalar exact_rho;

  // parameters
  Scalar x;
  Scalar t;

  //problem size
  int nx = 200;  // number of points
  int lx=10;     // length
  Scalar dx=Scalar(lx)/Scalar(nx);

  int nt = 100;  // number of points
  int lt = 10;     // length
  Scalar dt = Scalar(lt)/Scalar(nt);

  // initalize
  masa_init<Scalar>("euler-chemistry-test","euler_transient_1d");

  // initialize the default parameters
  masa_init_param<Scalar>();

  // check that all terms have been initialized
  masa_sanity_check<Scalar>();

  // evaluate MMS (1D+time)
  for(int i=0;i<nx;i++)
    for(int j=0;j<nt;j++)
      {
	t=j*dt;
	x=i*dx;
	
	// evalulate source terms
	ufield = masa_eval_source_rho_u  <Scalar>(x,t);
	efield = masa_eval_source_rho_e  <Scalar>(x,t);
	rhofield = masa_eval_source_rho  <Scalar>(x,t);

	exact_u = masa_eval_exact_u      <Scalar>(x,t);
	exact_p = masa_eval_exact_p      <Scalar>(x,t);
	exact_rho = masa_eval_exact_rho  <Scalar>(x,t);

      }
  
  return 0;

}

int main()
{
  int err=0;

  err += run_regression<double>();
  //err += run_regression<long double>();

  return err;
}

