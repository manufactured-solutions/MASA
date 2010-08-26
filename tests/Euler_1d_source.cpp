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
// $Author$
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <math.h>

double PI = acos(-1);

double SourceQ_rho (
  double x,
  double u_0,
  double u_x,
  double rho_0,
  double rho_x,
  double p_0,
  double p_x,
  double a_px,
  double a_rhox,
  double a_ux,
  double L)
{
  double Q_rho;
  Q_rho = rho_x * cos(a_rhox * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L)) * a_rhox * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L)) * a_ux * PI / L;
  return(Q_rho);
}

double SourceQ_u (
  double x,
  double u_0,
  double u_x,
  double rho_0,
  double rho_x,
  double p_0,
  double p_x,
  double a_px,
  double a_rhox,
  double a_ux,
  double L)
{
  double Q_u;
  Q_u = -sin(a_px * PI * x / L) * a_px * PI * p_x / L + rho_x * cos(a_rhox * PI * x / L) * pow(u_0 + u_x * sin(a_ux * PI * x / L), 0.2e1) * a_rhox * PI / L + 0.2e1 * u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L)) * (u_0 + u_x * sin(a_ux * PI * x / L)) * a_ux * PI / L;
  return(Q_u);
}

double SourceQ_e (
  double x,
  double u_0,
  double u_x,
  double rho_0,
  double rho_x,
  double p_0,
  double p_x,
  double a_px,
  double a_rhox,
  double a_ux,
  double Gamma,
  double mu,
  double L)
{
  double Q_e;
  Q_e = cos(a_rhox * PI * x / L) * rho_x * pow(u_0 + u_x * sin(a_ux * PI * x / L), 0.3e1) * a_rhox * PI / L / 0.2e1 + cos(a_ux * PI * x / L) * (p_0 + p_x * cos(a_px * PI * x / L)) * a_ux * PI * u_x * Gamma / L / (Gamma - 0.1e1) - Gamma * p_x * sin(a_px * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L)) * a_px * PI / L / (Gamma - 0.1e1) + 0.3e1 / 0.2e1 * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L)) * pow(u_0 + u_x * sin(a_ux * PI * x / L), 0.2e1) * a_ux * PI * u_x / L;
  return(Q_e);
}
