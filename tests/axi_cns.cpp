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
// $Author$
// $Id$
//
// axi_cns.cpp: program that tests axisymmetric navier stokes
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h> // for MASA_STRICT_REGRESSION
#include <cmath>
#include <masa.h>
#include <cmath>

#include <iostream>
#include <stdlib.h>

using namespace MASA;
using namespace std;

typedef double Real;

const Real pi = std::acos(Real(-1));
const Real threshold = 1.0e-15; // should be small enough to catch any obvious problems

Real nancheck(Real x)
{
  if(isnan(x))
    {
      cout << "MASA REGRESSION FAILURE:: nan found!\n";
      exit(1);
    }
  return 1;
}

Real anQ_p(Real r,Real z,Real p_0,Real p_1,Real rho_0,Real rho_1,Real u_1,Real w_0,Real w_1,Real a_pr,Real a_pz,Real a_rhor,Real a_rhoz,Real a_ur,Real a_uz,Real a_wr,Real a_wz,Real pi,Real L,Real Gamma,Real mu)
{
  Real exact_p = p_0 + p_1 * std::sin(a_pr * pi * r / L) * std::cos(a_pz * pi * z / L);
  return exact_p;
}
  
Real anQ_u (Real r,Real z,Real p_0,Real p_1,Real rho_0,Real rho_1,Real u_1,Real w_0,Real w_1,Real a_pr,Real a_pz,Real a_rhor,Real a_rhoz,Real a_ur,Real a_uz,Real a_wr,Real a_wz,Real pi,Real L,Real Gamma,Real mu)
{
  Real exact_u = u_1 * (std::cos(a_ur * pi * r / L) - 0.1e1) * std::sin(a_uz * pi * z / L);
  return exact_u;
} 
 
Real anQ_w (Real r,Real z,Real w_0,Real w_1,Real a_wr,Real a_wz,Real pi,Real L)
{
  Real exact_w = w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L);
  return exact_w;
}

Real anQ_rho (Real r,Real z,Real p_0,Real p_1,Real rho_0,Real rho_1,Real u_1,Real w_0,Real w_1,Real a_pr,Real a_pz,Real a_rhor,Real a_rhoz,Real a_ur,Real a_uz,Real a_wr,Real a_wz,Real pi,Real L,Real Gamma,Real mu)
{ 
  Real exact_rho = rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L);
  return exact_rho;
}

// ----------------------------------------
//   Source Terms
// ----------------------------------------

Real SourceQ_e(Real r,Real z,Real p_0,Real p_1,Real rho_0,Real rho_1,Real u_1,Real w_0,Real w_1,Real a_pr,Real a_pz,Real a_rhor,Real a_rhoz,Real a_ur,Real a_uz,Real a_wr,Real a_wz,Real pi,Real L,Real Gamma,Real mu,Real k, Real R)
{
  Real Q_e = -std::sin(a_ur * pi * r / L) * std::sin(a_uz * pi * z / L) * (0.3e1 * u_1 * u_1 * std::pow(cos(a_ur * pi * r / L) - 0.1e1, 0.2e1) * std::pow(sin(a_uz * pi * z / L), 0.2e1) + std::pow(w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L), 0.2e1)) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * a_ur * pi * u_1 / L / 0.2e1 + (u_1 * u_1 * std::pow(cos(a_ur * pi * r / L) - 0.1e1, 0.2e1) * std::pow(sin(a_uz * pi * z / L), 0.2e1) + std::pow(w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L), 0.2e1)) * std::cos(a_rhor * pi * r / L) * std::cos(a_rhoz * pi * z / L) * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * a_rhoz * pi * rho_1 / L / 0.2e1 + std::sin(a_pr * pi * r / L) * std::cos(a_pz * pi * z / L) * a_pr * a_pr * pi * pi * k * p_1 / R * std::pow(L, -0.2e1) / (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) + std::sin(a_pr * pi * r / L) * std::cos(a_pz * pi * z / L) * a_pz * a_pz * pi * pi * k * p_1 / R * std::pow(L, -0.2e1) / (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) - 0.2e1 * std::sin(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L) * std::cos(a_pr * pi * r / L) * std::cos(a_pz * pi * z / L) * a_rhor * a_pr * pi * pi * k * p_1 * rho_1 / R * std::pow(L, -0.2e1) * std::pow(rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L), -0.2e1) - 0.2e1 * std::cos(a_rhor * pi * r / L) * std::cos(a_rhoz * pi * z / L) * std::sin(a_pr * pi * r / L) * std::sin(a_pz * pi * z / L) * a_rhoz * a_pz * pi * pi * k * p_1 * rho_1 / R * std::pow(L, -0.2e1) * std::pow(rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L), -0.2e1) + 0.2e1 / 0.3e1 * std::sin(a_wr * pi * r / L) * std::cos(a_wz * pi * z / L) * u_1 * (std::cos(a_ur * pi * r / L) - 0.1e1) * std::sin(a_uz * pi * z / L) * a_wr * a_wz * pi * pi * mu * w_1 * std::pow(L, -0.2e1) + std::cos(a_pz * pi * z / L) * std::cos(a_pr * pi * r / L) * u_1 * (std::cos(a_ur * pi * r / L) - 0.1e1) * std::sin(a_uz * pi * z / L) * a_pr * pi * p_1 * Gamma / L / (Gamma - 0.1e1) - (u_1 * u_1 * std::pow(cos(a_ur * pi * r / L) - 0.1e1, 0.2e1) * std::pow(sin(a_uz * pi * z / L), 0.2e1) + std::pow(w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L), 0.2e1)) * std::sin(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L) * u_1 * (std::cos(a_ur * pi * r / L) - 0.1e1) * std::sin(a_uz * pi * z / L) * a_rhor * pi * rho_1 / L / 0.2e1 + std::cos(a_uz * pi * z / L) * std::pow(cos(a_ur * pi * r / L) - 0.1e1, 0.2e1) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * u_1 * u_1 * std::sin(a_uz * pi * z / L) * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * a_uz * pi / L - std::sin(a_ur * pi * r / L) * std::cos(a_uz * pi * z / L) * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * a_ur * a_uz * pi * pi * mu * u_1 * std::pow(L, -0.2e1) / 0.3e1 + std::cos(a_uz * pi * z / L) * (std::cos(a_ur * pi * r / L) - 0.1e1) * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * a_uz * pi * mu * u_1 / L / r / 0.3e1 - std::sin(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * u_1 * (std::cos(a_ur * pi * r / L) - 0.1e1) * std::sin(a_uz * pi * z / L) * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * a_wr * pi * w_1 / L - Gamma * std::sin(a_pr * pi * r / L) * std::sin(a_pz * pi * z / L) * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * a_pz * pi * p_1 / L / (Gamma - 0.1e1) - std::sin(a_ur * pi * r / L) * std::sin(a_uz * pi * z / L) * (p_0 + p_1 * std::sin(a_pr * pi * r / L) * std::cos(a_pz * pi * z / L)) * a_ur * pi * u_1 * Gamma / L / (Gamma - 0.1e1) - std::sin(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L) * (p_0 + p_1 * std::sin(a_pr * pi * r / L) * std::cos(a_pz * pi * z / L)) * a_rhor * pi * k * rho_1 / R / L * std::pow(rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L), -0.2e1) / r - 0.2e1 * std::pow(sin(a_rhoz * pi * z / L), 0.2e1) * (p_0 + p_1 * std::sin(a_pr * pi * r / L) * std::cos(a_pz * pi * z / L)) * std::pow(sin(a_rhor * pi * r / L), 0.2e1) * a_rhor * a_rhor * pi * pi * k * rho_1 * rho_1 / R * std::pow(L, -0.2e1) * std::pow(rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L), -0.3e1) - std::sin(a_rhoz * pi * z / L) * (p_0 + p_1 * std::sin(a_pr * pi * r / L) * std::cos(a_pz * pi * z / L)) * std::cos(a_rhor * pi * r / L) * a_rhor * a_rhor * pi * pi * k * rho_1 / R * std::pow(L, -0.2e1) * std::pow(rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L), -0.2e1) - 0.2e1 * std::pow(cos(a_rhor * pi * r / L), 0.2e1) * (p_0 + p_1 * std::sin(a_pr * pi * r / L) * std::cos(a_pz * pi * z / L)) * std::pow(cos(a_rhoz * pi * z / L), 0.2e1) * a_rhoz * a_rhoz * pi * pi * k * rho_1 * rho_1 / R * std::pow(L, -0.2e1) * std::pow(rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L), -0.3e1) - std::sin(a_rhoz * pi * z / L) * (p_0 + p_1 * std::sin(a_pr * pi * r / L) * std::cos(a_pz * pi * z / L)) * std::cos(a_rhor * pi * r / L) * a_rhoz * a_rhoz * pi * pi * k * rho_1 / R * std::pow(L, -0.2e1) * std::pow(rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L), -0.2e1) - a_uz * a_wr * pi * pi * mu * u_1 * w_1 * std::cos(a_uz * pi * z / L) * std::sin(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L) * (std::cos(a_ur * pi * r / L) - 0.1e1) * std::pow(L, -0.2e1) + 0.4e1 / 0.3e1 * a_ur * a_wz * pi * pi * mu * u_1 * w_1 * std::sin(a_ur * pi * r / L) * std::cos(a_wr * pi * r / L) * std::cos(a_wz * pi * z / L) * std::sin(a_uz * pi * z / L) * std::pow(L, -0.2e1) - std::cos(a_pz * pi * z / L) * std::cos(a_pr * pi * r / L) * a_pr * pi * k * p_1 / R / L / (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) / r - 0.4e1 / 0.3e1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L) * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * a_wz * a_wz * pi * pi * mu * w_1 * std::pow(L, -0.2e1) - 0.4e1 / 0.3e1 * std::pow(sin(a_uz * pi * z / L), 0.2e1) * u_1 * u_1 * (std::cos(a_ur * pi * r / L) - 0.1e1) * a_ur * a_ur * pi * pi * mu * std::pow(L, -0.2e1) + 0.4e1 / 0.3e1 * std::pow(cos(a_wr * pi * r / L), 0.2e1) * std::pow(cos(a_wz * pi * z / L), 0.2e1) * a_wz * a_wz * pi * pi * mu * w_1 * w_1 * std::pow(L, -0.2e1) + std::pow(cos(a_ur * pi * r / L) - 0.1e1, 0.2e1) * std::pow(cos(a_uz * pi * z / L), 0.2e1) * a_uz * a_uz * pi * pi * mu * u_1 * u_1 * std::pow(L, -0.2e1) + (p_0 + p_1 * std::sin(a_pr * pi * r / L) * std::cos(a_pz * pi * z / L)) * u_1 * (std::cos(a_ur * pi * r / L) - 0.1e1) * std::sin(a_uz * pi * z / L) * Gamma / (Gamma - 0.1e1) / r - 0.4e1 / 0.3e1 * u_1 * u_1 * std::pow(cos(a_ur * pi * r / L) - 0.1e1, 0.2e1) * std::pow(sin(a_uz * pi * z / L), 0.2e1) * a_ur * a_ur * pi * pi * mu * std::pow(L, -0.2e1) + 0.4e1 / 0.3e1 * std::pow(sin(a_uz * pi * z / L), 0.2e1) * std::pow(sin(a_ur * pi * r / L), 0.2e1) * a_ur * a_ur * pi * pi * mu * u_1 * u_1 * std::pow(L, -0.2e1) - u_1 * u_1 * std::pow(cos(a_ur * pi * r / L) - 0.1e1, 0.2e1) * std::pow(sin(a_uz * pi * z / L), 0.2e1) * a_uz * a_uz * pi * pi * mu * std::pow(L, -0.2e1) + (u_1 * u_1 * std::pow(cos(a_ur * pi * r / L) - 0.1e1, 0.2e1) * std::pow(sin(a_uz * pi * z / L), 0.2e1) + std::pow(w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L), 0.2e1)) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * u_1 * (std::cos(a_ur * pi * r / L) - 0.1e1) * std::sin(a_uz * pi * z / L) / r / 0.2e1 + a_wz * pi * w_1 * std::cos(a_wr * pi * r / L) * std::cos(a_wz * pi * z / L) * ((u_1 * u_1 * std::pow(cos(a_ur * pi * r / L) - 0.1e1, 0.2e1) * std::pow(sin(a_uz * pi * z / L), 0.2e1) + 0.3e1 * std::pow(w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L), 0.2e1)) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) / L / 0.2e1 - 0.4e1 / 0.3e1 * mu * u_1 * (std::cos(a_ur * pi * r / L) - 0.1e1) * std::sin(a_uz * pi * z / L) / L / r + Gamma * (p_0 + p_1 * std::sin(a_pr * pi * r / L) * std::cos(a_pz * pi * z / L)) / L / (Gamma - 0.1e1));
  return(Q_e);
}

Real SourceQ_u(Real r,Real z,Real p_0,Real p_1,Real rho_0,Real rho_1,Real u_1,Real w_0,Real w_1,Real a_pr,Real a_pz,Real a_rhor,Real a_rhoz,Real a_ur,Real a_uz,Real a_wr,Real a_wz,Real pi,Real L,Real Gamma,Real mu)
{
  Real Q_u = 0.4e1 / 0.3e1 * mu * u_1 * std::cos(a_ur * pi * r / L) * std::sin(a_uz * pi * z / L) * a_ur * a_ur * pi * pi * std::pow(L, -0.2e1) + mu * u_1 * (std::cos(a_ur * pi * r / L) - 0.1e1) * std::sin(a_uz * pi * z / L) * a_uz * a_uz * pi * pi * std::pow(L, -0.2e1) - 0.2e1 / 0.3e1 * mu * w_1 * std::sin(a_wr * pi * r / L) * std::cos(a_wz * pi * z / L) * a_wr * a_wz * pi * pi * std::pow(L, -0.2e1) + p_1 * std::cos(a_pr * pi * r / L) * std::cos(a_pz * pi * z / L) * a_pr * pi / L - u_1 * u_1 * std::pow(sin(a_uz * pi * z / L), 0.2e1) * std::pow(cos(a_ur * pi * r / L) - 0.1e1, 0.2e1) * rho_1 * std::sin(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L) * a_rhor * pi / L + u_1 * (std::cos(a_ur * pi * r / L) - 0.1e1) * rho_1 * std::cos(a_rhor * pi * r / L) * std::cos(a_rhoz * pi * z / L) * std::sin(a_uz * pi * z / L) * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * a_rhoz * pi / L + u_1 * (std::cos(a_ur * pi * r / L) - 0.1e1) * std::cos(a_uz * pi * z / L) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * a_uz * pi / L + 0.2e1 / 0.3e1 * std::sin(a_ur * pi * r / L) * std::sin(a_uz * pi * z / L) * a_ur * pi * u_1 * mu / L / r - 0.2e1 * std::sin(a_ur * pi * r / L) * std::pow(sin(a_uz * pi * z / L), 0.2e1) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * u_1 * u_1 * (std::cos(a_ur * pi * r / L) - 0.1e1) * a_ur * pi / L + 0.2e1 / 0.3e1 * std::cos(a_wr * pi * r / L) * std::cos(a_wz * pi * z / L) * a_wz * pi * w_1 * mu / L / r + std::cos(a_wr * pi * r / L) * std::cos(a_wz * pi * z / L) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * u_1 * (std::cos(a_ur * pi * r / L) - 0.1e1) * std::sin(a_uz * pi * z / L) * a_wz * pi * w_1 / L + u_1 * u_1 * std::pow(sin(a_uz * pi * z / L), 0.2e1) * std::pow(cos(a_ur * pi * r / L) - 0.1e1, 0.2e1) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) / r;
  return(Q_u);  
}

Real SourceQ_w(Real r,Real z,Real p_0,Real p_1,Real rho_0,Real rho_1,Real u_1,Real w_0,Real w_1,Real a_pr,Real a_pz,Real a_rhor,Real a_rhoz,Real a_ur,Real a_uz,Real a_wr,Real a_wz,Real pi,Real L,Real Gamma,Real mu)
{
  Real Q_w = std::sin(a_ur * pi * r / L) * std::cos(a_uz * pi * z / L) * a_ur * a_uz * pi * pi * mu * u_1 * std::pow(L, -0.2e1) / 0.3e1 + 0.4e1 / 0.3e1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L) * a_wz * a_wz * pi * pi * mu * w_1 * std::pow(L, -0.2e1) - std::sin(a_pr * pi * r / L) * std::sin(a_pz * pi * z / L) * a_pz * pi * p_1 / L - std::sin(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L) * u_1 * (std::cos(a_ur * pi * r / L) - 0.1e1) * std::sin(a_uz * pi * z / L) * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * a_rhor * pi * rho_1 / L + std::cos(a_rhor * pi * r / L) * std::cos(a_rhoz * pi * z / L) * std::pow(w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L), 0.2e1) * a_rhoz * pi * rho_1 / L - std::sin(a_uz * pi * z / L) * std::sin(a_ur * pi * r / L) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * a_ur * pi * u_1 / L - std::sin(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * u_1 * (std::cos(a_ur * pi * r / L) - 0.1e1) * std::sin(a_uz * pi * z / L) * a_wr * pi * w_1 / L + 0.2e1 * std::cos(a_wz * pi * z / L) * std::cos(a_wr * pi * r / L) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * a_wz * pi * w_1 / L - (std::cos(a_ur * pi * r / L) - 0.1e1) * std::cos(a_uz * pi * z / L) * a_uz * pi * mu * u_1 / L / r / 0.3e1 + std::sin(a_uz * pi * z / L) * u_1 * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * (std::cos(a_ur * pi * r / L) - 0.1e1) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) / r;
  return(Q_w);
}

Real SourceQ_rho(Real r,Real z,Real p_0,Real p_1,Real rho_0,Real rho_1,Real u_1,Real w_0,Real w_1,Real a_pr,Real a_pz,Real a_rhor,Real a_rhoz,Real a_ur,Real a_uz,Real a_wr,Real a_wz,Real pi,Real L,Real Gamma,Real mu)
{
  Real Q_rho = -u_1 * std::sin(a_uz * pi * z / L) * rho_1 * std::sin(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L) * (std::cos(a_ur * pi * r / L) - 0.1e1) * a_rhor * pi / L + rho_1 * std::cos(a_rhor * pi * r / L) * std::cos(a_rhoz * pi * z / L) * (w_0 + w_1 * std::cos(a_wr * pi * r / L) * std::sin(a_wz * pi * z / L)) * a_rhoz * pi / L - u_1 * std::sin(a_uz * pi * z / L) * std::sin(a_ur * pi * r / L) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * a_ur * pi / L + w_1 * std::cos(a_wr * pi * r / L) * std::cos(a_wz * pi * z / L) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * a_wz * pi / L + u_1 * std::sin(a_uz * pi * z / L) * (rho_0 + rho_1 * std::cos(a_rhor * pi * r / L) * std::sin(a_rhoz * pi * z / L)) * (std::cos(a_ur * pi * r / L) - 0.1e1) / r;
  return(Q_rho);
}

int main()
{  
  //variables
  Real R;  
  Real p_0;
  Real p_1;
  Real rho_0;
  Real rho_1;
  Real u_1;
  Real w_0;
  Real w_1;
  Real a_pr;
  Real a_pz;
  Real a_rhor;
  Real a_rhoz;
  Real a_ur;
  Real a_uz;
  Real a_wr;
  Real a_wz;
  Real L;
  Real mu;
  Real Gamma;    
  Real k;
  
  // parameters
  Real r;
  Real z;

  // solutions
  Real ufield,ufield2,ufield3;
  Real wfield,wfield2,wfield3;
  Real efield,efield2,efield3;
  Real rho,rho2,rho3;

  Real exact_u,exact_u2,exact_u3;
  Real exact_w,exact_w2,exact_w3;
  Real exact_p,exact_p2,exact_p3;
  Real exact_rho,exact_rho2,exact_rho3;

  // initalize
  int nx = 115;  // number of points
  int ny = 68;  
  int lx=3;     // length
  int ly=1; 
  
  Real dx=Real(lx)/Real(nx);
  Real dy=Real(ly)/Real(ny);

  masa_init<Real>("axisymmetric_navierstokes_compressible","axisymmetric_navierstokes_compressible");

  // set params
  masa_init_param<Real>();
  
  // get vars
  R      = masa_get_param<Real>("R");
  p_0    = masa_get_param<Real>("p_0");
  p_1    = masa_get_param<Real>("p_1");
  rho_0  = masa_get_param<Real>("rho_0");
  rho_1  = masa_get_param<Real>("rho_1");
  u_1    = masa_get_param<Real>("u_1");
  w_0    = masa_get_param<Real>("w_0");
  w_1    = masa_get_param<Real>("w_1");
  a_pr   = masa_get_param<Real>("a_pr");
  a_pz   = masa_get_param<Real>("a_pz");
  a_rhor = masa_get_param<Real>("a_rhor");
  a_rhoz = masa_get_param<Real>("a_rhoz");
  a_ur   = masa_get_param<Real>("a_ur");
  a_uz   = masa_get_param<Real>("a_uz");
  a_wr   = masa_get_param<Real>("a_wr");
  a_wz   = masa_get_param<Real>("a_wz");
  L      = masa_get_param<Real>("L");
  Gamma  = masa_get_param<Real>("Gamma");
  mu     = masa_get_param<Real>("mu");
  k      =masa_get_param<Real>("k");

  // check that all terms have been initialized
  int err = masa_sanity_check<Real>();
  if(err != 0)
    {
      cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }
  
  // evaluate source terms (2D)
  for(int i=1;i<nx;i++)
    for(int j=1;j<ny;j++)    
      {
	r=i*dx;
	z=j*dy;
	
	//evalulate source terms
	ufield = masa_eval_source_rho_u<Real>  (r,z);
	wfield = masa_eval_source_rho_w<Real>  (r,z);
	efield = masa_eval_source_rho_e<Real>  (r,z);
	rho     = masa_eval_source_rho<Real>(r,z);

	//evaluate analytical terms
	exact_u = masa_eval_exact_u<Real>        (r,z);
	exact_w = masa_eval_exact_w<Real>        (r,z);
	exact_p = masa_eval_exact_p<Real>        (r,z);
	exact_rho = masa_eval_exact_rho<Real>    (r,z);
	  
	// check against maple
	ufield2 = SourceQ_u   (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma, mu);
	wfield2 = SourceQ_w   (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma, mu);
	rho2    = SourceQ_rho (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma, mu);
	efield2 = SourceQ_e   (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma, mu,k,R);
	
	exact_u2   = anQ_u   (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma, mu);
	exact_w2   = anQ_w   (r, z, w_0, w_1, a_wr, a_wz, pi, L);
	exact_rho2 = anQ_rho (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma, mu);
	exact_p2   = anQ_p   (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma, mu);

	// test the result is roughly zero
	// choose between abs and rel error
#ifdef MASA_STRICT_REGRESSION

	ufield3 = std::abs(ufield-ufield2);
	wfield3 = std::abs(wfield-wfield2);
	efield3 = std::abs(efield-efield2);
	rho3    = std::abs(rho-rho2);
	
	exact_u3   = std::abs(exact_u-exact_u2);
	exact_w3  = std::abs(exact_w-exact_w2);
	exact_rho3 = std::abs(exact_rho-exact_rho2);
	exact_p3   = std::abs(exact_p-exact_p2);

#else
	ufield3 = std::abs(ufield-ufield2)/std::abs(ufield2);
	wfield3 = std::abs(wfield-wfield2)/std::abs(wfield2);
	efield3 = std::abs(efield-efield2)/std::abs(efield2);
	rho3    = std::abs(rho-rho2)/std::abs(rho2);
	
	exact_u3   = std::abs(exact_u-exact_u2)/std::abs(exact_u2);
	exact_w3   = std::abs(exact_w-exact_w2)/std::abs(exact_w2);
	exact_rho3 = std::abs(exact_rho-exact_rho2)/std::abs(exact_rho2);
	exact_p3   = std::abs(exact_p-exact_p2)/std::abs(exact_p2);
#endif

	nancheck(ufield3);
	nancheck(wfield3);
	nancheck(efield3);
	nancheck(rho3);

	nancheck(exact_u3);
	nancheck(exact_w3);
	nancheck(exact_rho3);
	nancheck(exact_p3);

	if(ufield3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Navier-Stokes\n";
	    cout << "U Field Source Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << ufield3 << endl;
	    cout << "Source term is:                   " << ufield2 << endl;
	    cout << "MASA term is:                     " << ufield << endl;
	    cout << r << " " << z << endl;
	    exit(1);
	  }

	if(exact_u3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Navier-Stokes\n";
	    cout << "U Field Analytical Term\n";
	    cout << "Exceeded Threshold by: " << exact_u << endl;
	    cout.precision(16);
	    cout << r << " " << z << endl;
	    exit(1);
	  }

	if(wfield3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Navier-Stokes\n";
	    cout << "W Field Source Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by:            " << wfield3 << endl;
	    cout << "Source term is:                   " << wfield2 << endl;
	    cout << "MASA term is:                     " << wfield << endl;
	    cout << r << " " << z << endl;
	    exit(1);
	  }

	// this guy is broken
	if(exact_w3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Navier-Stokes\n";
	    cout << "W Field Analytical Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << exact_w3 << endl;
	    cout << "Source term is:        " << exact_w2 << endl;
	    cout << "MASA term is:          " << exact_w << endl;
	    cout << r << " " << z << endl;
	    exit(1);
	  }

	if(efield3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Navier-Stokes\n";
	    cout << "Energy Source Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << efield3 << endl;
	    cout << "Source term is:        " << efield2 << endl;
	    cout << "MASA term is:          " << efield << endl;
	    cout << r << " " << z << endl;
	    exit(1);
	  }

	if(exact_p3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Navier-Stokes\n";
	    cout << "P Field Analytical Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << exact_p << endl;
	    cout << r << " " << z << endl;
	    exit(1);
	  }

	if(rho3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Navier-Stokes\n";
	    cout.precision(16);
	    cout << "RHO Source Term\n";
	    cout << "Exceeded Threshold by: " << rho << endl;
	    cout << r << " " << z << endl;
	    exit(1);
	  }

	if(exact_rho3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Navier-Stokes\n";
	    cout.precision(16);
	    cout << "RHO Analytical Term\n";
	    cout << "Exceeded Threshold by: " << exact_rho << endl;
	    cout << r << " " << z << endl;
	    exit(1);
	  }

      } // done iterating

  // tests passed
  return 0;
}
