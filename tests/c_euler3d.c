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
// c_euler3d.c :program that tests masa against known source term
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h> // for MASA_STRICT_REGRESSION
#include <masa.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

double anQ_p (double x,double y,double z,double p_0,double p_x,double p_y,double p_z,double a_px,double a_py,double a_pz,double L)
{
  double pi = acos(-1);
  double p_an = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L);
  return p_an;
}
 
double anQ_u (double x,double y,double z,double u_0,double u_x,double u_y,double u_z,double a_ux,double a_uy,double a_uz,double L)
{
  double pi = acos(-1); 
  double u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L);  
  return u_an;
} 
 
double anQ_v (double x,double y,double z,double v_0,double v_x,double v_y,double v_z,double a_vx,double a_vy,double a_vz,double L)
{
  double pi = acos(-1);
  double v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L);
  return v_an;
}

double anQ_w (double x,double y,double z,double w_0,double w_x,double w_y,double w_z,double a_wx,double a_wy,double a_wz,double L)
{
  double pi = acos(-1);
  double w_an = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L);
  return w_an;
}

double anQ_rho (double x,double y,double z,double rho_0,double rho_x,double rho_y,double rho_z,double a_rhox,double a_rhoy,double a_rhoz,double L)
{ 
  double pi = acos(-1);
  double rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L);
  return rho_an;
}

double SourceQ_e ( // 40
  double x,
  double y,
  double z,
  double u_0,
  double u_x,
  double u_y,
  double u_z,
  double v_0,
  double v_x,
  double v_y,
  double v_z,
  double w_0,
  double w_x,
  double w_y,
  double w_z,
  double rho_0,
  double rho_x,
  double rho_y,
  double rho_z,
  double p_0,
  double p_x,
  double p_y,
  double p_z,
  double a_px,
  double a_py,
  double a_pz,
  double a_rhox,
  double a_rhoy,
  double a_rhoz,
  double a_ux,
  double a_uy,
  double a_uz,
  double a_vx,
  double a_vy,
  double a_vz,
  double a_wx,
  double a_wy,
  double a_wz,
  double mu,
  double Gamma,
  double L)
{
  double PI = acos(-1);
  double Q_e;
  Q_e = -Gamma * p_x * sin(a_px * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_px * PI / L / (Gamma - 0.1e1) + Gamma * p_y * cos(a_py * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_py * PI / L / (Gamma - 0.1e1) - Gamma * p_z * sin(a_pz * PI * z / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_pz * PI / L / (Gamma - 0.1e1) + (pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) + pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1)) * a_rhox * PI * rho_x * cos(a_rhox * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) / L / 0.2e1 - (pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) + pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1)) * a_rhoy * PI * rho_y * sin(a_rhoy * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) / L / 0.2e1 + (pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) + pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1)) * a_rhoz * PI * rho_z * cos(a_rhoz * PI * z / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) / L / 0.2e1 + a_ux * PI * u_x * cos(a_ux * PI * x / L) * ((0.3e1 * rho_x * sin(a_rhox * PI * x / L) + 0.3e1 * rho_y * cos(a_rhoy * PI * y / L) + 0.3e1 * rho_z * sin(a_rhoz * PI * z / L) + 0.3e1 * rho_0) * pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) + (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L), 0.2e1) + (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1)) / L / 0.2e1 + Gamma * a_ux * PI * u_x * cos(a_ux * PI * x / L) * (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L) + p_0) / L / (Gamma - 0.1e1) - (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * PI * a_uy * u_y * sin(a_uy * PI * y / L) / L - (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * PI * a_uz * u_z * sin(a_uz * PI * z / L) / L - (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * PI * a_vx * v_x * sin(a_vx * PI * x / L) / L + a_vy * PI * v_y * cos(a_vy * PI * y / L) * ((rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) + (0.3e1 * rho_x * sin(a_rhox * PI * x / L) + 0.3e1 * rho_y * cos(a_rhoy * PI * y / L) + 0.3e1 * rho_z * sin(a_rhoz * PI * z / L) + 0.3e1 * rho_0) * pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L), 0.2e1) + (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1)) / L / 0.2e1 + Gamma * a_vy * PI * v_y * cos(a_vy * PI * y / L) * (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L) + p_0) / L / (Gamma - 0.1e1) + (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * PI * a_vz * v_z * cos(a_vz * PI * z / L) / L + (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * PI * a_wx * w_x * cos(a_wx * PI * x / L) / L + (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * PI * a_wy * w_y * cos(a_wy * PI * y / L) / L - a_wz * PI * w_z * sin(a_wz * PI * z / L) * ((rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) + (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L), 0.2e1) + (0.3e1 * rho_x * sin(a_rhox * PI * x / L) + 0.3e1 * rho_y * cos(a_rhoy * PI * y / L) + 0.3e1 * rho_z * sin(a_rhoz * PI * z / L) + 0.3e1 * rho_0) * pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1)) / L / 0.2e1 - Gamma * a_wz * PI * w_z * sin(a_wz * PI * z / L) * (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L) + p_0) / L / (Gamma - 0.1e1);
  return(Q_e);
}

double SourceQ_u ( // 38 variables
  double x,
  double y,
  double z,
  double u_0,
  double u_x,
  double u_y,
  double u_z,
  double v_0,
  double v_x,
  double v_y,
  double v_z,
  double w_0,
  double w_x,
  double w_y,
  double w_z,
  double rho_0,
  double rho_x,
  double rho_y,
  double rho_z,
  double p_0,
  double p_x,
  double p_y,
  double p_z,
  double a_px,
  double a_py,
  double a_pz,
  double a_rhox,
  double a_rhoy,
  double a_rhoz,
  double a_ux,
  double a_uy,
  double a_uz,
  double a_vx,
  double a_vy,
  double a_vz,
  double a_wx,
  double a_wy,
  double a_wz,
  double L)
{
  double PI = acos(-1);
  double Q_u;
  Q_u = -p_x * sin(a_px * PI * x / L) * a_px * PI / L + rho_x * cos(a_rhox * PI * x / L) * pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_rhoy * PI / L + rho_z * cos(a_rhoz * PI * z / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_rhoz * PI / L + 0.2e1 * u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_ux * PI / L - u_y * sin(a_uy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_uy * PI / L - u_z * sin(a_uz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_uz * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_vy * PI / L - w_z * sin(a_wz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_wz * PI / L;
  return(Q_u);
}

double SourceQ_v (
  double x,
  double y,
  double z,
  double u_0,
  double u_x,
  double u_y,
  double u_z,
  double v_0,
  double v_x,
  double v_y,
  double v_z,
  double w_0,
  double w_x,
  double w_y,
  double w_z,
  double rho_0,
  double rho_x,
  double rho_y,
  double rho_z,
  double p_0,
  double p_x,
  double p_y,
  double p_z,
  double a_px,
  double a_py,
  double a_pz,
  double a_rhox,
  double a_rhoy,
  double a_rhoz,
  double a_ux,
  double a_uy,
  double a_uz,
  double a_vx,
  double a_vy,
  double a_vz,
  double a_wx,
  double a_wy,
  double a_wz,
  double L)
{
  double PI = acos(-1);
  double Q_v;
  Q_v = p_y * cos(a_py * PI * y / L) * a_py * PI / L + rho_x * cos(a_rhox * PI * x / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L), 0.2e1) * a_rhoy * PI / L + rho_z * cos(a_rhoz * PI * z / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_rhoz * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_ux * PI / L - v_x * sin(a_vx * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_vx * PI / L + 0.2e1 * v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_vy * PI / L + v_z * cos(a_vz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_vz * PI / L - w_z * sin(a_wz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_wz * PI / L;
  return(Q_v);
}

double SourceQ_w (
  double x,
  double y,
  double z,
  double u_0,
  double u_x,
  double u_y,
  double u_z,
  double v_0,
  double v_x,
  double v_y,
  double v_z,
  double w_0,
  double w_x,
  double w_y,
  double w_z,
  double rho_0,
  double rho_x,
  double rho_y,
  double rho_z,
  double p_0,
  double p_x,
  double p_y,
  double p_z,
  double a_px,
  double a_py,
  double a_pz,
  double a_rhox,
  double a_rhoy,
  double a_rhoz,
  double a_ux,
  double a_uy,
  double a_uz,
  double a_vx,
  double a_vy,
  double a_vz,
  double a_wx,
  double a_wy,
  double a_wz,
  double L)
{
  double PI = acos(-1);
  double Q_w;
  Q_w = -p_z * sin(a_pz * PI * z / L) * a_pz * PI / L + rho_x * cos(a_rhox * PI * x / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_rhoy * PI / L + rho_z * cos(a_rhoz * PI * z / L) * pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1) * a_rhoz * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_ux * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_vy * PI / L + w_x * cos(a_wx * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_wx * PI / L + w_y * cos(a_wy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_wy * PI / L - 0.2e1 * w_z * sin(a_wz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_wz * PI / L;
  return(Q_w);
}

double SourceQ_rho(
  double x,
  double y,
  double z,
  double u_0,
  double u_x,
  double u_y,
  double u_z,
  double v_0,
  double v_x,
  double v_y,
  double v_z,
  double w_0,
  double w_x,
  double w_y,
  double w_z,
  double rho_0,
  double rho_x,
  double rho_y,
  double rho_z,
  double p_0,
  double p_x,
  double p_y,
  double p_z,
  double a_px,
  double a_py,
  double a_pz,
  double a_rhox,
  double a_rhoy,
  double a_rhoz,
  double a_ux,
  double a_uy,
  double a_uz,
  double a_vx,
  double a_vy,
  double a_vz,
  double a_wx,
  double a_wy,
  double a_wz,
  double mu,
  double L)
{
  double PI = acos(-1);
  double Q_rho;
  Q_rho = rho_x * cos(a_rhox * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_rhoy * PI / L + rho_z * cos(a_rhoz * PI * z / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_rhoz * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * a_ux * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * a_vy * PI / L - w_z * sin(a_wz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * a_wz * PI / L;
  return(Q_rho);
}

int main()
{
  //variables
  double u_0;
  double u_x;
  double u_y;
  double u_z;
  double v_0;
  double v_x;
  double v_y;
  double v_z;
  double w_0;
  double w_x;
  double w_y;
  double w_z;
  double rho_0;
  double rho_x;
  double rho_y;
  double rho_z;
  double p_0;
  double p_x;
  double p_y;
  double p_z;
  double a_px;
  double a_py;
  double a_pz;
  double a_rhox;
  double a_rhoy;
  double a_rhoz;
  double a_ux;
  double a_uy;
  double a_uz;
  double a_vx;
  double a_vy;
  double a_vz;
  double a_wx;
  double a_wy;
  double a_wz;
  double mu;
  double Gamma;
  double L;    

  // parameters
  double x;
  double y;
  double z;
  int i,j,k;

  // solutions
  double ufield,ufield2,ufield3;
  double vfield,vfield2,vfield3;
  double wfield,wfield2,wfield3;
  double efield,efield2,efield3;
  double rho,rho2,rho3;

  double u_an,u_an2,u_an3;
  double v_an,v_an2,v_an3;
  double w_an,w_an2,w_an3;
  double p_an,p_an2,p_an3;
  double rho_an,rho_an2,rho_an3;

  // initalize
  int nx = 43;             // number of points
  int ny = 12;  
  int nz = 23;  
  int lx=1;                // length
  int ly=2; 
  int lz=3;   
  double dx=(double)(lx)/(double)(nx);
  double dy=(double)(ly)/(double)(ny);
  double dz=(double)(lz)/(double)(nz);

  cmasa_init("euler-test","euler_3d");

  // set params
  cmasa_init_param();

  // get vars for comparison
  u_0 = cmasa_get_param("u_0");
  u_x = cmasa_get_param("u_x");
  u_y = cmasa_get_param("u_y");
  u_z = cmasa_get_param("u_z");

  v_0 = cmasa_get_param("v_0");
  v_x = cmasa_get_param("v_x");
  v_y = cmasa_get_param("v_y");
  v_z = cmasa_get_param("v_z");

  w_0 = cmasa_get_param("w_0");
  w_x = cmasa_get_param("w_x");
  w_y = cmasa_get_param("w_y");
  w_z = cmasa_get_param("w_z");

  rho_0 = cmasa_get_param("rho_0");
  rho_x = cmasa_get_param("rho_x");
  rho_y = cmasa_get_param("rho_y");
  rho_z = cmasa_get_param("rho_z");

  p_0 = cmasa_get_param("p_0");
  p_x = cmasa_get_param("p_x");
  p_y = cmasa_get_param("p_y");
  p_z = cmasa_get_param("p_z");

  a_px = cmasa_get_param("a_px");
  a_py = cmasa_get_param("a_py");
  a_pz = cmasa_get_param("a_pz");

  a_rhox = cmasa_get_param("a_rhox");
  a_rhoy = cmasa_get_param("a_rhoy");
  a_rhoz = cmasa_get_param("a_rhoz");

  a_ux = cmasa_get_param("a_ux");
  a_uy = cmasa_get_param("a_uy");
  a_uz = cmasa_get_param("a_uz");

  a_vx = cmasa_get_param("a_vx");
  a_vy = cmasa_get_param("a_vy");
  a_vz = cmasa_get_param("a_vz");

  a_wx = cmasa_get_param("a_wx");
  a_wy = cmasa_get_param("a_wy");
  a_wz = cmasa_get_param("a_wz");

  Gamma = cmasa_get_param("Gamma");
  mu    = cmasa_get_param("mu");
  L     = cmasa_get_param("L");

  // check all vars initialized
  cmasa_sanity_check();

  // evaluate source terms (3D)
  for(i=0;i<nx;i++)
    for(j=0;j<ny;j++)    
      for(k=0;k<nz;k++)
	{
	  x=i*dx;
	  y=j*dy;
	  z=k*dz;

	  //evalulate source terms
	  ufield = cmasa_eval_3d_u_source  (x,y,z);
	  vfield = cmasa_eval_3d_v_source  (x,y,z);
	  wfield = cmasa_eval_3d_w_source  (x,y,z);
	  efield = cmasa_eval_3d_e_source  (x,y,z);
	  rho    = cmasa_eval_3d_rho_source(x,y,z);
	  
	  //evaluate analytical terms
	  u_an   = cmasa_eval_3d_u_an      (x,y,z);
	  v_an   = cmasa_eval_3d_v_an      (x,y,z);
	  w_an   = cmasa_eval_3d_w_an      (x,y,z);
	  p_an   = cmasa_eval_3d_p_an      (x,y,z);
	  rho_an = cmasa_eval_3d_rho_an     (x,y,z);

	  // check against maple output
	  ufield2   = SourceQ_u  (x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,L);
	  vfield2   = SourceQ_v  (x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,L);
	  wfield2   = SourceQ_w  (x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,L);
	  rho2      = SourceQ_rho(x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,mu,L);
	  efield2   = SourceQ_e  (x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,mu,Gamma,L);
  
	  u_an2     = anQ_u   (x,y,z,u_0,u_x,u_y,u_z,a_ux,a_uy,a_uz,L);
	  v_an2     = anQ_v   (x,y,z,v_0,v_x,v_y,v_z,a_vx,a_vy,a_vz,L);
	  w_an2     = anQ_w   (x,y,z,w_0,w_x,w_y,w_z,a_wx,a_wy,a_wz,L);
	  rho_an2   = anQ_rho (x,y,z,rho_0,rho_x,rho_y,rho_z,a_rhox,a_rhoy,a_rhoz,L);
	  p_an2     = anQ_p   (x,y,z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,L);

	  // test the result is roughly zero
	  // choose between abs and rel error
#ifdef MASA_STRICT_REGRESSION

	  ufield3 = fabs(ufield-ufield2);
	  vfield3 = fabs(vfield-vfield2);
	  wfield3 = fabs(wfield-wfield2);
	  efield3 = fabs(efield-efield2);
	  rho3    = fabs(rho-rho2);

	  u_an3   = fabs(u_an-u_an2);
	  v_an3   = fabs(v_an-v_an2);
	  w_an3   = fabs(w_an-w_an2);
	  rho_an3 = fabs(rho_an-rho_an2);
	  p_an3   = fabs(p_an-p_an2);

#else

	  ufield3 = fabs(ufield-ufield2)/fabs(ufield2);
	  vfield3 = fabs(vfield-vfield2)/fabs(vfield2);
	  wfield3 = fabs(wfield-wfield2)/fabs(wfield2);
	  efield3 = fabs(efield-efield2)/fabs(efield2);
	  rho3    = fabs(rho-rho2)/fabs(rho2);

	  u_an3   = fabs(u_an-u_an2)/fabs(u_an2);
	  v_an3   = fabs(v_an-v_an2)/fabs(v_an2);
	  w_an3   = fabs(w_an-w_an2)/fabs(w_an2);
	  rho_an3 = fabs(rho_an-rho_an2)/fabs(rho_an2);
	  p_an3   = fabs(p_an-p_an2)/fabs(p_an2);

#endif

	  if(ufield3 > threshold)
	    {
	      printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-3d\n");
	      printf("U Field Source Term\n");
	      exit(1);
	    }
	  
	  if(u_an3 > threshold)
	    {
	      printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-3d\n");
	      printf("U Field Analytical Term\n");
	      exit(1);
	    }

	  if(vfield3 > threshold)
	    {
	      printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-3d\n");
	      printf("V Field Source Term\n");
	      exit(1);
	    }

	  if(v_an3 > threshold)
	    {
	      printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-3d\n");
	      printf("V Field Analytical Term\n");
	      exit(1);
	    }

	  if(wfield3 > threshold)
	    {
	      printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-3d\n");
	      printf("W Field Source Term\n");
	      printf("Threshold Exceeded: %g\n",wfield);
	      exit(1);
	    }

	  if(w_an3 > threshold)
	    {
	      printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-3d\n");
	      printf("W Field Analytical Term\n");
	      printf("Threshold Exceeded: %g\n\n",w_an);
	      exit(1);
	    }

	  if(efield3 > threshold)
	    {
	      printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-3d\n");
	      printf("E Field Source Term\n");
	      printf("Threshold Exceeded: %g\n",efield3);
	      printf("CMASA:              %5.16f\n",efield);
	      printf("Maple:              %5.16f\n",efield2);
	      printf("x,y:                %g %g\n\n",x,y);	   
	      exit(1);
	    }

	  if(p_an3 > threshold)
	    {
	      printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-3d\n");
	      printf("P Field Analytical Term\n");
	      exit(1);
	    }

	  if(rho3 > threshold)
	    {
	      printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-3d\n");
	      printf("RHO Field Source Term\n");
	      exit(1);
	    }

	  if(rho_an3 > threshold)
	    {
	      printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-3d\n");
	      printf("RHO Field Analytical Term\n");
	      exit(1);
	    }
	  
	}// done iterating

  // tests passed
  return 0;
}
