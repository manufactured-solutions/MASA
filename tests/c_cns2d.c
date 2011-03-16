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
//
// program that tests navier-stokes-2d against known source term
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h> // for MASA_STRICT_REGRESSION
#include <masa.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef MASA_STRICT_REGRESSION
const double threshold = 1.0e-15;
#else
const double threshold = 2.0e-15; // should be small enough to catch any obvious problems
#endif

double anQ_p (double x,double y,double p_0,double p_x,double p_y,double a_px,double a_py,double L)
{
  const double pi = acos(-1);
  double exact_p = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L);
  return exact_p;
}
  
double anQ_u (double x,double y,double u_0,double u_x,double u_y,double a_ux,double a_uy,double L)
{
  const double pi = acos(-1);
 
  double exact_u = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  return exact_u;
} 
 
double anQ_v (double x,double y,double v_0,double v_x,double v_y,double a_vx,double a_vy,double L)
{
  const double pi = acos(-1);

  double exact_v = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);
  return exact_v;
}

double anQ_rho (double x,double y,double rho_0,double rho_x,double rho_y,double a_rhox,double a_rhoy,double L)
{ 
  const double pi = acos(-1);
 
  double exact_rho = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  return exact_rho;
}


double SourceQ_u (
  double x,
  double y,
  double u_0,
  double u_x,
  double u_y,
  double v_0,
  double v_x,
  double v_y,
  double rho_0,
  double rho_x,
  double rho_y,
  double p_0,
  double p_x,
  double p_y,
  double a_px,
  double a_py,
  double a_rhox,
  double a_rhoy,
  double a_ux,
  double a_uy,
  double a_vx,
  double a_vy,
  double mu,
  double L,
  double R,
  double k)
{
  const double PI = acos(-1);
  double Q_u;
  Q_u = (double)(0.4e1) / (double)0.3e1 * mu * u_x * sin(a_ux * PI * x / L) * a_ux * a_ux * PI * PI * pow(L, -0.2e1) + mu * u_y * cos(a_uy * PI * y / L) * a_uy * a_uy * PI * PI * pow(L, -0.2e1) - p_x * sin(a_px * PI * x / L) * a_px * PI / L + rho_x * cos(a_rhox * PI * x / L) * pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L), 0.2e1) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_rhoy * PI / L + 0.2e1 * u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_ux * PI / L - u_y * sin(a_uy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_uy * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_vy * PI / L;
  return(Q_u);
}

double SourceQ_v (
  double x,
  double y,
  double u_0,
  double u_x,
  double u_y,
  double v_0,
  double v_x,
  double v_y,
  double rho_0,
  double rho_x,
  double rho_y,
  double p_0,
  double p_x,
  double p_y,
  double a_px,
  double a_py,
  double a_rhox,
  double a_rhoy,
  double a_ux,
  double a_uy,
  double a_vx,
  double a_vy,
  double mu,
  double L,
  double R,
  double k)
{
  const double PI = acos(-1);
  double Q_v;
  Q_v = mu * v_x * cos(a_vx * PI * x / L) * a_vx * a_vx * PI * PI * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * mu * v_y * sin(a_vy * PI * y / L) * a_vy * a_vy * PI * PI * pow(L, -0.2e1) + p_y * cos(a_py * PI * y / L) * a_py * PI / L + rho_x * cos(a_rhox * PI * x / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L), 0.2e1) * a_rhoy * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_ux * PI / L - v_x * sin(a_vx * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_vx * PI / L + 0.2e1 * v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_vy * PI / L;
  return(Q_v);
}

#include <math.h>

double SourceQ_e (
  double x,
  double y,
  double u_0,
  double u_x,
  double u_y,
  double v_0,
  double v_x,
  double v_y,
  double rho_0,
  double rho_x,
  double rho_y,
  double p_0,
  double p_x,
  double p_y,
  double a_px,
  double a_py,
  double a_rhox,
  double a_rhoy,
  double a_ux,
  double a_uy,
  double a_vx,
  double a_vy,
  double Gamma,
  double mu,
  double L,
  double R,
  double k)
{
  const double pi = acos(-1);
  double Q_e;
  Q_e = -(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * (pow(u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0, (double)2) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0, (double)2)) * rho_y * sin(a_rhoy * pi * y / L) * a_rhoy * pi / L / (double)2 + (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * (pow(u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0, (double)2) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0, (double)2)) * rho_x * cos(a_rhox * pi * x / L) * a_rhox * pi / L / (double)2 + (double)4 / (double)3 * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * mu * v_y * sin(a_vy * pi * y / L) * a_vy * a_vy * pi * pi * pow(L, -(double)2) - (double)4 / (double)3 * mu * v_y * v_y * pow(cos(a_vy * pi * y / L), (double)2) * a_vy * a_vy * pi * pi * pow(L, -(double)2) - mu * v_x * v_x * pow(sin(a_vx * pi * x / L), (double)2) * a_vx * a_vx * pi * pi * pow(L, -(double)2) - (double)4 / (double)3 * mu * u_x * u_x * pow(cos(a_ux * pi * x / L), (double)2) * a_ux * a_ux * pi * pi * pow(L, -(double)2) - mu * u_y * u_y * pow(sin(a_uy * pi * y / L), (double)2) * a_uy * a_uy * pi * pi * pow(L, -(double)2) + (Gamma * (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_0) / (Gamma - (double)1) + (pow(u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0, (double)2) / (double)2 + (double)3 / (double)2 * pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0, (double)2)) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0)) * v_y * cos(a_vy * pi * y / L) * a_vy * pi / L + (Gamma * (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_0) / (Gamma - (double)1) + ((double)3 / (double)2 * pow(u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0, (double)2) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0, (double)2) / (double)2) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0)) * u_x * cos(a_ux * pi * x / L) * a_ux * pi / L + (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * mu * v_x * cos(a_vx * pi * x / L) * a_vx * a_vx * pi * pi * pow(L, -(double)2) + (double)4 / (double)3 * (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * mu * u_x * sin(a_ux * pi * x / L) * a_ux * a_ux * pi * pi * pow(L, -(double)2) + (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * mu * u_y * cos(a_uy * pi * y / L) * a_uy * a_uy * pi * pi * pow(L, -(double)2) - (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * u_y * sin(a_uy * pi * y / L) * a_uy * pi / L - (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_0) * rho_x * k * sin(a_rhox * pi * x / L) * a_rhox * a_rhox * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -(double)2) * pow(L, -(double)2) / R - ((double)2 * p_x * cos(a_px * pi * x / L) + (double)2 * p_y * sin(a_py * pi * y / L) + (double)2 * p_0) * rho_x * rho_x * k * pow(cos(a_rhox * pi * x / L), (double)2) * a_rhox * a_rhox * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -(double)3) * pow(L, -(double)2) / R - (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_0) * rho_y * k * cos(a_rhoy * pi * y / L) * a_rhoy * a_rhoy * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -(double)2) * pow(L, -(double)2) / R - ((double)2 * p_x * cos(a_px * pi * x / L) + (double)2 * p_y * sin(a_py * pi * y / L) + (double)2 * p_0) * rho_y * rho_y * k * pow(sin(a_rhoy * pi * y / L), (double)2) * a_rhoy * a_rhoy * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -(double)3) * pow(L, -(double)2) / R + (double)4 / (double)3 * mu * u_x * v_y * cos(a_ux * pi * x / L) * cos(a_vy * pi * y / L) * a_ux * a_vy * pi * pi * pow(L, -(double)2) - (double)2 * mu * u_y * v_x * sin(a_uy * pi * y / L) * sin(a_vx * pi * x / L) * a_uy * a_vx * pi * pi * pow(L, -(double)2) - (double)2 * k * p_x * rho_x * cos(a_rhox * pi * x / L) * sin(a_px * pi * x / L) * a_px * a_rhox * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -(double)2) * pow(L, -(double)2) / R - (double)2 * k * p_y * rho_y * cos(a_py * pi * y / L) * sin(a_rhoy * pi * y / L) * a_py * a_rhoy * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -(double)2) * pow(L, -(double)2) / R - (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * v_x * sin(a_vx * pi * x / L) * a_vx * pi / L - Gamma * (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * p_x * sin(a_px * pi * x / L) * a_px * pi / (Gamma - (double)1) / L + Gamma * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * p_y * cos(a_py * pi * y / L) * a_py * pi / (Gamma - (double)1) / L + k * p_x * cos(a_px * pi * x / L) * a_px * a_px * pi * pi / (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * pow(L, -(double)2) / R + k * p_y * sin(a_py * pi * y / L) * a_py * a_py * pi * pi / (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * pow(L, -(double)2) / R;


  //Q_e = -(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * (pow(u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1) + pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0, 0.2e1)) * rho_y * sin(a_rhoy * PI * y / L) * a_rhoy * PI / L / 0.2e1 + (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * (pow(u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1) + pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0, 0.2e1)) * rho_x * cos(a_rhox * PI * x / L) * a_rhox * PI / L / 0.2e1 + 0.4e1 / 0.3e1 * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * mu * v_y * sin(a_vy * PI * y / L) * a_vy * a_vy * PI * PI * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * mu * v_y * v_y * pow(cos(a_vy * PI * y / L), 0.2e1) * a_vy * a_vy * PI * PI * pow(L, -0.2e1) - mu * v_x * v_x * pow(sin(a_vx * PI * x / L), 0.2e1) * a_vx * a_vx * PI * PI * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * mu * u_x * u_x * pow(cos(a_ux * PI * x / L), 0.2e1) * a_ux * a_ux * PI * PI * pow(L, -0.2e1) - mu * u_y * u_y * pow(sin(a_uy * PI * y / L), 0.2e1) * a_uy * a_uy * PI * PI * pow(L, -0.2e1) + (Gamma * (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_0) / (Gamma - 0.1e1) + (pow(u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1) / 0.2e1 + 0.3e1 / 0.2e1 * pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0, 0.2e1)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0)) * v_y * cos(a_vy * PI * y / L) * a_vy * PI / L + (Gamma * (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_0) / (Gamma - 0.1e1) + (0.3e1 / 0.2e1 * pow(u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1) + pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0, 0.2e1) / 0.2e1) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0)) * u_x * cos(a_ux * PI * x / L) * a_ux * PI / L + (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * mu * v_x * cos(a_vx * PI * x / L) * a_vx * a_vx * PI * PI * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * mu * u_x * sin(a_ux * PI * x / L) * a_ux * a_ux * PI * PI * pow(L, -0.2e1) + (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * mu * u_y * cos(a_uy * PI * y / L) * a_uy * a_uy * PI * PI * pow(L, -0.2e1) - (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * u_y * sin(a_uy * PI * y / L) * a_uy * PI / L - (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_0) * rho_x * k * sin(a_rhox * PI * x / L) * a_rhox * a_rhox * PI * PI * pow(rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0, -0.2e1) * pow(L, -0.2e1) / R - (0.2e1 * p_x * cos(a_px * PI * x / L) + 0.2e1 * p_y * sin(a_py * PI * y / L) + 0.2e1 * p_0) * rho_x * rho_x * k * pow(cos(a_rhox * PI * x / L), 0.2e1) * a_rhox * a_rhox * PI * PI * pow(rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0, -0.3e1) * pow(L, -0.2e1) / R - (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_0) * rho_y * k * cos(a_rhoy * PI * y / L) * a_rhoy * a_rhoy * PI * PI * pow(rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0, -0.2e1) * pow(L, -0.2e1) / R - (0.2e1 * p_x * cos(a_px * PI * x / L) + 0.2e1 * p_y * sin(a_py * PI * y / L) + 0.2e1 * p_0) * rho_y * rho_y * k * pow(sin(a_rhoy * PI * y / L), 0.2e1) * a_rhoy * a_rhoy * PI * PI * pow(rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0, -0.3e1) * pow(L, -0.2e1) / R + 0.4e1 / 0.3e1 * mu * u_x * v_y * cos(a_ux * PI * x / L) * cos(a_vy * PI * y / L) * a_ux * a_vy * PI * PI * pow(L, -0.2e1) - 0.2e1 * mu * u_y * v_x * sin(a_uy * PI * y / L) * sin(a_vx * PI * x / L) * a_uy * a_vx * PI * PI * pow(L, -0.2e1) - 0.2e1 * k * p_x * rho_x * cos(a_rhox * PI * x / L) * sin(a_px * PI * x / L) * a_px * a_rhox * PI * PI * pow(rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0, -0.2e1) * pow(L, -0.2e1) / R - 0.2e1 * k * p_y * rho_y * cos(a_py * PI * y / L) * sin(a_rhoy * PI * y / L) * a_py * a_rhoy * PI * PI * pow(rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0, -0.2e1) * pow(L, -0.2e1) / R - (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * v_x * sin(a_vx * PI * x / L) * a_vx * PI / L - Gamma * (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * p_x * sin(a_px * PI * x / L) * a_px * PI / (Gamma - 0.1e1) / L + Gamma * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * p_y * cos(a_py * PI * y / L) * a_py * PI / (Gamma - 0.1e1) / L + k * p_x * cos(a_px * PI * x / L) * a_px * a_px * PI * PI / (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * pow(L, -0.2e1) / R + k * p_y * sin(a_py * PI * y / L) * a_py * a_py * PI * PI / (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * pow(L, -0.2e1) / R;
  return(Q_e);
}

double SourceQ_rho (
  double x,
  double y,
  double u_0,
  double u_x,
  double u_y,
  double v_0,
  double v_x,
  double v_y,
  double rho_0,
  double rho_x,
  double rho_y,
  double p_0,
  double p_x,
  double p_y,
  double a_px,
  double a_py,
  double a_rhox,
  double a_rhoy,
  double a_ux,
  double a_uy,
  double a_vx,
  double a_vy,
  double mu,
  double L,
  double R,
  double k)
{
  const double PI = acos(-1);
  double Q_rho;
  Q_rho = (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * a_rhox * PI * rho_x * cos(a_rhox * PI * x / L) / L - (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * a_rhoy * PI * rho_y * sin(a_rhoy * PI * y / L) / L + (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * a_ux * PI * u_x * cos(a_ux * PI * x / L) / L + (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * a_vy * PI * v_y * cos(a_vy * PI * y / L) / L;
  return(Q_rho);
}

int main()
{
  //variables
  double u_0;
  double u_x;
  double u_y;
  double v_0;
  double v_x;
  double v_y;
  double rho_0;
  double rho_x;
  double rho_y;
  double p_0;
  double p_x;
  double p_y;
  double a_px;
  double a_py;
  double a_rhox;
  double a_rhoy;
  double a_ux;
  double a_uy;
  double a_vx;
  double a_vy;
  double Gamma;
  double mu;
  double L;
  double R;
  double k;

  // parameters
  double x;
  double y;
  double z;
  int i,j;

  // solutions
  double ufield,ufield2,ufield3;
  double vfield,vfield2,vfield3;
  double wfield,wfield2,wfield3;
  double efield,efield2,efield3;
  double rho,rho2,rho3;

  double exact_u,exact_u2,exact_u3;
  double exact_v,exact_v2,exact_v3;
  double exact_w,exact_w2,exact_w3;
  double exact_p,exact_p2,exact_p3;
  double exact_rho,exact_rho2,exact_rho3;

  // initalize
  int nx = 10;  // number of points
  int ny = 8;  
  int lx = 2;     // length
  int ly = 1; 
    
  double dx=(double)(lx)/(double)(nx);
  double dy=(double)(ly)/(double)(ny);

  masa_init("navier-stokes-test","navierstokes_2d_compressible");

  // set params
  masa_init_param();
  
  // get vars for comparison
  u_0 = masa_get_param("u_0");
  u_x = masa_get_param("u_x");
  u_y = masa_get_param("u_y");

  v_0 = masa_get_param("v_0");
  v_x = masa_get_param("v_x");
  v_y = masa_get_param("v_y");

  rho_0 = masa_get_param("rho_0");
  rho_x = masa_get_param("rho_x");
  rho_y = masa_get_param("rho_y");

  p_0 = masa_get_param("p_0");
  p_x = masa_get_param("p_x");
  p_y = masa_get_param("p_y");

  a_px = masa_get_param("a_px");
  a_py = masa_get_param("a_py");

  a_rhox = masa_get_param("a_rhox");
  a_rhoy = masa_get_param("a_rhoy");

  a_ux = masa_get_param("a_ux");
  a_uy = masa_get_param("a_uy");

  a_vx = masa_get_param("a_vx");
  a_vy = masa_get_param("a_vy");

  Gamma = masa_get_param("Gamma");
  mu    = masa_get_param("mu");
  L     = masa_get_param("L");

  R = masa_get_param("R");
  k = masa_get_param("k");

  // check that all terms have been initialized
  masa_sanity_check();

  // evaluate source terms (2D)
  for(i=0;i<nx;i++)
    for(j=0;j<ny;j++)    
      {
	x=i*dx;
	y=j*dy;

	//evalulate source terms
	ufield = masa_eval_2d_source_rho_u  (x,y);
	vfield = masa_eval_2d_source_rho_v  (x,y);
	efield = masa_eval_2d_source_rho_e  (x,y);
	rho    = masa_eval_2d_source_rho(x,y);
	
	//evaluate analytical terms
	exact_u   = masa_eval_2d_exact_u      (x,y);
	exact_v   = masa_eval_2d_exact_v      (x,y);
	exact_p   = masa_eval_2d_exact_p      (x,y);
	exact_rho = masa_eval_2d_exact_rho    (x,y);

	// check against maple
	ufield2   = SourceQ_u  (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,mu,L,R,k);
	vfield2   = SourceQ_v  (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,mu,L,R,k);
	rho2      = SourceQ_rho(x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,mu,L,R,k);
	efield2   = SourceQ_e  (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,Gamma,mu,L,R,k);
	
	exact_u2   = anQ_u   (x,y,u_0,u_x,u_y,a_ux,a_uy,L);
	exact_v2   = anQ_v   (x,y,v_0,v_x,v_y,a_vx,a_vy,L);
	exact_rho2 = anQ_rho (x,y,rho_0,rho_x,rho_y,a_rhox,a_rhoy,L);
	exact_p2   = anQ_p   (x,y,p_0,p_x,p_y,a_px,a_py,L);

	// test the result is roughly zero
	// choose between abs and rel error
#ifdef MASA_STRICT_REGRESSION
	ufield3 = fabs(ufield-ufield2);
	vfield3 = fabs(vfield-vfield2);
	efield3 = fabs(efield-efield2);
	rho3    = fabs(rho-rho2);

	exact_u3   = fabs(exact_u-exact_u2);
	exact_v3   = fabs(exact_v-exact_v2);
	exact_rho3 = fabs(exact_rho-exact_rho2);
	exact_p3   = fabs(exact_p-exact_p2);
#else
	ufield3 = fabs(ufield-ufield2)/fabs(ufield2);
	vfield3 = fabs(vfield-vfield2)/fabs(vfield2);
	efield3 = fabs(efield-efield2)/fabs(efield2);
	rho3    = fabs(rho-rho2)/fabs(rho2);

	exact_u3   = fabs(exact_u-exact_u2)/fabs(exact_u2);
	exact_v3   = fabs(exact_v-exact_v2)/fabs(exact_v2);
	exact_rho3 = fabs(exact_rho-exact_rho2)/fabs(exact_rho2);
	exact_p3   = fabs(exact_p-exact_p2)/fabs(exact_p2);
#endif

	if(ufield3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Navier-Stokes 2d\n");
	    printf("U Field Source Term\n");
	    exit(1);
	  }

	if(exact_u3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Navier-Stokes 2d\n");
	    printf("U Field Analytical Term\n");
	    exit(1);
	  }

	if(vfield3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Navier-Stokes 2d\n");
	    printf("V Field Source Term\n");
	    printf("(Relative) Threshold Exceeded: %g\n",vfield3);
	    printf("MASA:              %5.16f\n",vfield);
	    printf("Maple:              %5.16f\n",vfield2);
	    printf("x,y:                %g %g\n",x,y);

	    exit(1);
	  }

	if(exact_v3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Navier-Stokes 2d\n");
	    printf("V Field Analytical Term\n");
	    exit(1);
	  }

	if(efield3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Navier-Stokes 2d\n");
	    printf("E Field Source Term\n");
	    printf("(Relative) Threshold Exceeded: %g\n",efield3);
	    printf("MASA:              %5.16f\n",efield);
	    printf("Maple:              %5.16f\n",efield2);
	    printf("x,y:                %g %g\n",x,y);
	    exit(1);
	  }

	if(exact_p3 > threshold)
	  {	    
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Navier-Stokes 2d\n");
	    printf("P Field Analytical Term\n");
	    exit(1);
	  }

	if(rho3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Navier-Stokes 2d\n");
	    printf("RHO Field Source Term\n");
	    exit(1);
	  }

	if(exact_rho3 > threshold)
	  {	    
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Navier-Stokes 2d\n");
	    printf("RHO Analytical Term\n");
	    exit(1);
	  }

      } // done iterating

  // tests passed
  return 0;
}
