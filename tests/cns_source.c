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
// $Id: euler-source.c 12693 2010-08-26 03:35:34Z nick $
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <math.h>

double eval_2d_exact_p(double x,double y,double p_0,double p_x,double p_y,double a_px,double a_py,double L)
{

  const double pi = acos(-1);
 
  double exact_p = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L);
  return exact_p;
}
  
double eval_2d_exact_u(double x,double y,double u_0,double u_x,double u_y,double a_ux,double a_uy,double L)
{
  const double pi = acos(-1);

  double exact_u = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  return exact_u;
} 
 
double eval_2d_exact_v(double x,double y,double v_0,double v_x,double v_y,double a_vx,double a_vy,double L)
{
  const double pi = acos(-1);

  double exact_v = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);
  return exact_v;
}

double eval_2d_exact_rho(double x,double y,double rho_0,double rho_x,double rho_y,double a_rhox,double a_rhoy,double L)
{ 
  const double pi = acos(-1);

  double exact_rho = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  return exact_rho;
}

double eval_2d_source_e(
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
  Q_e = -(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * (pow(u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0, 0.2e1) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0, 0.2e1)) * rho_y * sin(a_rhoy * pi * y / L) * a_rhoy * pi / L / 0.2e1 + (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * (pow(u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0, 0.2e1) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0, 0.2e1)) * rho_x * cos(a_rhox * pi * x / L) * a_rhox * pi / L / 0.2e1 + 0.4e1 / 0.3e1 * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * mu * v_y * sin(a_vy * pi * y / L) * a_vy * a_vy * pi * pi * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * mu * v_y * v_y * pow(cos(a_vy * pi * y / L), 0.2e1) * a_vy * a_vy * pi * pi * pow(L, -0.2e1) - mu * v_x * v_x * pow(sin(a_vx * pi * x / L), 0.2e1) * a_vx * a_vx * pi * pi * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * mu * u_x * u_x * pow(cos(a_ux * pi * x / L), 0.2e1) * a_ux * a_ux * pi * pi * pow(L, -0.2e1) - mu * u_y * u_y * pow(sin(a_uy * pi * y / L), 0.2e1) * a_uy * a_uy * pi * pi * pow(L, -0.2e1) + (Gamma * (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_0) / (Gamma - 0.1e1) + (pow(u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0, 0.2e1) / 0.2e1 + 0.3e1 / 0.2e1 * pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0, 0.2e1)) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0)) * v_y * cos(a_vy * pi * y / L) * a_vy * pi / L + (Gamma * (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_0) / (Gamma - 0.1e1) + (0.3e1 / 0.2e1 * pow(u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0, 0.2e1) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0, 0.2e1) / 0.2e1) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0)) * u_x * cos(a_ux * pi * x / L) * a_ux * pi / L + (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * mu * v_x * cos(a_vx * pi * x / L) * a_vx * a_vx * pi * pi * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * mu * u_x * sin(a_ux * pi * x / L) * a_ux * a_ux * pi * pi * pow(L, -0.2e1) + (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * mu * u_y * cos(a_uy * pi * y / L) * a_uy * a_uy * pi * pi * pow(L, -0.2e1) - (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * u_y * sin(a_uy * pi * y / L) * a_uy * pi / L - (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_0) * rho_x * k * sin(a_rhox * pi * x / L) * a_rhox * a_rhox * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -0.2e1) * pow(L, -0.2e1) / R - (0.2e1 * p_x * cos(a_px * pi * x / L) + 0.2e1 * p_y * sin(a_py * pi * y / L) + 0.2e1 * p_0) * rho_x * rho_x * k * pow(cos(a_rhox * pi * x / L), 0.2e1) * a_rhox * a_rhox * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -0.3e1) * pow(L, -0.2e1) / R - (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_0) * rho_y * k * cos(a_rhoy * pi * y / L) * a_rhoy * a_rhoy * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -0.2e1) * pow(L, -0.2e1) / R - (0.2e1 * p_x * cos(a_px * pi * x / L) + 0.2e1 * p_y * sin(a_py * pi * y / L) + 0.2e1 * p_0) * rho_y * rho_y * k * pow(sin(a_rhoy * pi * y / L), 0.2e1) * a_rhoy * a_rhoy * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -0.3e1) * pow(L, -0.2e1) / R + 0.4e1 / 0.3e1 * mu * u_x * v_y * cos(a_ux * pi * x / L) * cos(a_vy * pi * y / L) * a_ux * a_vy * pi * pi * pow(L, -0.2e1) - 0.2e1 * mu * u_y * v_x * sin(a_uy * pi * y / L) * sin(a_vx * pi * x / L) * a_uy * a_vx * pi * pi * pow(L, -0.2e1) - 0.2e1 * k * p_x * rho_x * cos(a_rhox * pi * x / L) * sin(a_px * pi * x / L) * a_px * a_rhox * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -0.2e1) * pow(L, -0.2e1) / R - 0.2e1 * k * p_y * rho_y * cos(a_py * pi * y / L) * sin(a_rhoy * pi * y / L) * a_py * a_rhoy * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -0.2e1) * pow(L, -0.2e1) / R - (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * v_x * sin(a_vx * pi * x / L) * a_vx * pi / L - Gamma * (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * p_x * sin(a_px * pi * x / L) * a_px * pi / (Gamma - 0.1e1) / L + Gamma * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * p_y * cos(a_py * pi * y / L) * a_py * pi / (Gamma - 0.1e1) / L + k * p_x * cos(a_px * pi * x / L) * a_px * a_px * pi * pi / (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * pow(L, -0.2e1) / R + k * p_y * sin(a_py * pi * y / L) * a_py * a_py * pi * pi / (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * pow(L, -0.2e1) / R;
  return(Q_e);
}

double eval_2d_source_u(
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
  const double pi = acos(-1);

  double Q_u;
  Q_u = 0.4e1 / 0.3e1 * mu * u_x * sin(a_ux * pi * x / L) * a_ux * a_ux * pi * pi * pow(L, -0.2e1) + mu * u_y * cos(a_uy * pi * y / L) * a_uy * a_uy * pi * pi * pow(L, -0.2e1) - p_x * sin(a_px * pi * x / L) * a_px * pi / L + rho_x * cos(a_rhox * pi * x / L) * pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L), 0.2e1) * a_rhox * pi / L - rho_y * sin(a_rhoy * pi * y / L) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_rhoy * pi / L + 0.2e1 * u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_ux * pi / L - u_y * sin(a_uy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * a_uy * pi / L + v_y * cos(a_vy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_vy * pi / L;
  return(Q_u);
}

double eval_2d_source_v(
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
			 double L)
{
  const double pi = acos(-1);

  double Q_v;
  Q_v = mu * v_x * cos(a_vx * pi * x / L) * a_vx * a_vx * pi * pi * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * mu * v_y * sin(a_vy * pi * y / L) * a_vy * a_vy * pi * pi * pow(L, -0.2e1) + p_y * cos(a_py * pi * y / L) * a_py * pi / L + rho_x * cos(a_rhox * pi * x / L) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_rhox * pi / L - rho_y * sin(a_rhoy * pi * y / L) * pow(v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L), 0.2e1) * a_rhoy * pi / L + u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * a_ux * pi / L - v_x * sin(a_vx * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_vx * pi / L + 0.2e1 * v_y * cos(a_vy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * a_vy * pi / L;
  return(Q_v);
}

double eval_2d_source_rho(
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
  const double pi = acos(-1);

  double Q_rho;
  Q_rho = (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * a_rhox * pi * rho_x * cos(a_rhox * pi * x / L) / L - (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * a_rhoy * pi * rho_y * sin(a_rhoy * pi * y / L) / L + (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * a_ux * pi * u_x * cos(a_ux * pi * x / L) / L + (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * a_vy * pi * v_y * cos(a_vy * pi * y / L) / L;
  return(Q_rho);
}

//--------------------------------------------------------------------------
// 
// 
//    3D problem
// 
//
//--------------------------------------------------------------------------

double eval_3d_exact_p(double x,double y,double z,double p_0,double p_x,double p_y,double p_z,double a_px,double a_py,double a_pz,double L)
{
  const double pi = acos(-1);
  double exact_p = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L);
  return exact_p;
}

double eval_3d_exact_u(double x,double y,double z,double u_0,double u_x,double u_y,double u_z,double a_ux,double a_uy,double a_uz,double L)
{
  const double pi = acos(-1);
  double exact_u = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L);  
  return exact_u;
} 

double eval_3d_exact_v(double x,double y,double z,double v_0,double v_x,double v_y,double v_z,double a_vx,double a_vy,double a_vz,double L)
{
  const double pi = acos(-1);
  double exact_v = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L);
  return exact_v;
}
double eval_3d_exact_w(double x,double y,double z,double w_0,double w_x,double w_y,double w_z,double a_wx,double a_wy,double a_wz,double L)
{
  const double pi = acos(-1);
  double exact_w = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L);
  return exact_w;
}
double eval_3d_exact_rho(double x,double y,double z,double rho_0,double rho_x,double rho_y,double rho_z,double a_rhox,double a_rhoy,double a_rhoz,double L)
{ 
  const double pi = acos(-1);
  double exact_rho = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L);
  return exact_rho;
}

double eval_3d_source_e(  
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
  double L,
  double R,
  double k)
{    
  const double pi = acos(-1);
  // onkars generated routines
  double Q_e = cos(a_rhox * pi * x / L) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * (pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), 0.2e1) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0, 0.2e1)) * rho_x * a_rhox * pi / L / 0.2e1 - sin(a_rhoy * pi * y / L) * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0) * (pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), 0.2e1) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0, 0.2e1)) * rho_y * a_rhoy * pi / L / 0.2e1 + cos(a_rhoz * pi * z / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), 0.2e1) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0, 0.2e1)) * rho_z * a_rhoz * pi / L / 0.2e1 + ((pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), 0.2e1) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0, 0.2e1) + 0.3e1 * pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), 0.2e1)) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) / L / 0.2e1 + Gamma * (p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L)) / L / (Gamma - 0.1e1)) * u_x * cos(a_ux * pi * x / L) * a_ux * pi + ((pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), 0.2e1) + 0.3e1 * pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0, 0.2e1)) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) / L / 0.2e1 + Gamma * (p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L)) / L / (Gamma - 0.1e1)) * v_y * cos(a_vy * pi * y / L) * a_vy * pi + (-(0.3e1 * pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), 0.2e1) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0, 0.2e1) + pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), 0.2e1)) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) / L / 0.2e1 - Gamma * (p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L)) / L / (Gamma - 0.1e1)) * w_z * sin(a_wz * pi * z / L) * a_wz * pi + 0.4e1 / 0.3e1 * (-pow(cos(a_ux * pi * x / L), 0.2e1) * u_x + sin(a_ux * pi * x / L) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L))) * mu * u_x * a_ux * a_ux * pi * pi * pow(L, -0.2e1) + (-pow(sin(a_uy * pi * y / L), 0.2e1) * u_y + cos(a_uy * pi * y / L) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L))) * mu * u_y * a_uy * a_uy * pi * pi * pow(L, -0.2e1) + (-pow(sin(a_uz * pi * z / L), 0.2e1) * u_z + cos(a_uz * pi * z / L) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L))) * mu * u_z * a_uz * a_uz * pi * pi * pow(L, -0.2e1) - (pow(sin(a_vx * pi * x / L), 0.2e1) * v_x - cos(a_vx * pi * x / L) * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0)) * mu * v_x * a_vx * a_vx * pi * pi * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * (pow(cos(a_vy * pi * y / L), 0.2e1) * v_y - sin(a_vy * pi * y / L) * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0)) * mu * v_y * a_vy * a_vy * pi * pi * pow(L, -0.2e1) - (pow(cos(a_vz * pi * z / L), 0.2e1) * v_z - sin(a_vz * pi * z / L) * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0)) * mu * v_z * a_vz * a_vz * pi * pi * pow(L, -0.2e1) + (-pow(cos(a_wx * pi * x / L), 0.2e1) * w_x + sin(a_wx * pi * x / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L))) * mu * w_x * a_wx * a_wx * pi * pi * pow(L, -0.2e1) + (-pow(cos(a_wy * pi * y / L), 0.2e1) * w_y + sin(a_wy * pi * y / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L))) * mu * w_y * a_wy * a_wy * pi * pi * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * (-pow(sin(a_wz * pi * z / L), 0.2e1) * w_z + cos(a_wz * pi * z / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L))) * mu * w_z * a_wz * a_wz * pi * pi * pow(L, -0.2e1) + sin(a_py * pi * y / L) * k * p_y * a_py * a_py * pi * pi * pow(L, -0.2e1) / R / (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) - 0.2e1 * cos(a_rhox * pi * x / L) * rho_x * sin(a_px * pi * x / L) * k * p_x * a_px * a_rhox * pi * pi * pow(L, -0.2e1) / R * pow(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L), -0.2e1) - 0.2e1 * sin(a_rhoy * pi * y / L) * rho_y * cos(a_py * pi * y / L) * k * p_y * a_py * a_rhoy * pi * pi * pow(L, -0.2e1) / R * pow(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L), -0.2e1) - (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * u_y * sin(a_uy * pi * y / L) * a_uy * pi / L + (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0) * v_z * cos(a_vz * pi * z / L) * a_vz * pi / L + cos(a_px * pi * x / L) * k * p_x * a_px * a_px * pi * pi * pow(L, -0.2e1) / R / (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) + cos(a_pz * pi * z / L) * k * p_z * a_pz * a_pz * pi * pi * pow(L, -0.2e1) / R / (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) + (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * w_x * cos(a_wx * pi * x / L) * a_wx * pi / L - (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * v_x * sin(a_vx * pi * x / L) * a_vx * pi / L - (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * u_z * sin(a_uz * pi * z / L) * a_uz * pi / L + (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0) * w_y * cos(a_wy * pi * y / L) * a_wy * pi / L - 0.2e1 * cos(a_rhoz * pi * z / L) * rho_z * sin(a_pz * pi * z / L) * k * p_z * a_pz * a_rhoz * pi * pi * pow(L, -0.2e1) / R * pow(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L), -0.2e1) - (0.2e1 * pow(cos(a_rhox * pi * x / L), 0.2e1) * rho_x + sin(a_rhox * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L))) * (p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L)) * k * rho_x * a_rhox * a_rhox * pi * pi * pow(L, -0.2e1) / R * pow(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L), -0.3e1) - (0.2e1 * pow(sin(a_rhoy * pi * y / L), 0.2e1) * rho_y + cos(a_rhoy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L))) * (p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L)) * k * rho_y * a_rhoy * a_rhoy * pi * pi * pow(L, -0.2e1) / R * pow(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L), -0.3e1) - (0.2e1 * pow(cos(a_rhoz * pi * z / L), 0.2e1) * rho_z + sin(a_rhoz * pi * z / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L))) * (p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L)) * k * rho_z * a_rhoz * a_rhoz * pi * pi * pow(L, -0.2e1) / R * pow(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L), -0.3e1) + 0.4e1 / 0.3e1 * mu * u_x * v_y * cos(a_ux * pi * x / L) * cos(a_vy * pi * y / L) * a_ux * a_vy * pi * pi * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * mu * u_x * w_z * cos(a_ux * pi * x / L) * sin(a_wz * pi * z / L) * a_ux * a_wz * pi * pi * pow(L, -0.2e1) - 0.2e1 * mu * u_y * v_x * sin(a_uy * pi * y / L) * sin(a_vx * pi * x / L) * a_uy * a_vx * pi * pi * pow(L, -0.2e1) + 0.2e1 * mu * u_z * w_x * cos(a_wx * pi * x / L) * sin(a_uz * pi * z / L) * a_uz * a_wx * pi * pi * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * mu * v_y * w_z * cos(a_vy * pi * y / L) * sin(a_wz * pi * z / L) * a_vy * a_wz * pi * pi * pow(L, -0.2e1) - 0.2e1 * mu * v_z * w_y * cos(a_vz * pi * z / L) * cos(a_wy * pi * y / L) * a_vz * a_wy * pi * pi * pow(L, -0.2e1) - Gamma * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * sin(a_px * pi * x / L) * p_x * a_px * pi / L / (Gamma - 0.1e1) + Gamma * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0) * cos(a_py * pi * y / L) * p_y * a_py * pi / L / (Gamma - 0.1e1) - Gamma * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * sin(a_pz * pi * z / L) * p_z * a_pz * pi / L / (Gamma - 0.1e1);  
  return(Q_e);
}

double eval_3d_source_u(  
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
  double L,
  double R,
  double k)
{
  const double pi = acos(-1);
  double Q_u;
  Q_u = 0.4e1 / 0.3e1 * mu * u_x * sin(a_ux * pi * x / L) * a_ux * a_ux * pi * pi * pow(L, -0.2e1) + mu * u_y * cos(a_uy * pi * y / L) * a_uy * a_uy * pi * pi * pow(L, -0.2e1) + mu * u_z * cos(a_uz * pi * z / L) * a_uz * a_uz * pi * pi * pow(L, -0.2e1) - p_x * sin(a_px * pi * x / L) * a_px * pi / L + rho_x * cos(a_rhox * pi * x / L) * pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), 0.2e1) * a_rhox * pi / L - rho_y * sin(a_rhoy * pi * y / L) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_rhoy * pi / L + rho_z * cos(a_rhoz * pi * z / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_rhoz * pi / L + 0.2e1 * u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_ux * pi / L - u_y * sin(a_uy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_uy * pi / L - u_z * sin(a_uz * pi * z / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * a_uz * pi / L + v_y * cos(a_vy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_vy * pi / L - w_z * sin(a_wz * pi * z / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_wz * pi / L;
  return(Q_u);
}

double eval_3d_source_v(  
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
  double L,
  double R,
  double k)
{
  const double pi = acos(-1);
  double Q_v;
  Q_v = mu * v_x * cos(a_vx * pi * x / L) * a_vx * a_vx * pi * pi * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * mu * v_y * sin(a_vy * pi * y / L) * a_vy * a_vy * pi * pi * pow(L, -0.2e1) + mu * v_z * sin(a_vz * pi * z / L) * a_vz * a_vz * pi * pi * pow(L, -0.2e1) + p_y * cos(a_py * pi * y / L) * a_py * pi / L + rho_x * cos(a_rhox * pi * x / L) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_rhox * pi / L - rho_y * sin(a_rhoy * pi * y / L) * pow(v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L), 0.2e1) * a_rhoy * pi / L + rho_z * cos(a_rhoz * pi * z / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_rhoz * pi / L + u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_ux * pi / L - v_x * sin(a_vx * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_vx * pi / L + 0.2e1 * v_y * cos(a_vy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_vy * pi / L + v_z * cos(a_vz * pi * z / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * a_vz * pi / L - w_z * sin(a_wz * pi * z / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_wz * pi / L;
  return(Q_v);
}

double eval_3d_source_w(  double x,
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
  double L,
  double R,
  double k)
{
  const double pi = acos(-1);
  double Q_w;
  Q_w = mu * w_x * sin(a_wx * pi * x / L) * a_wx * a_wx * pi * pi * pow(L, -0.2e1) + mu * w_y * sin(a_wy * pi * y / L) * a_wy * a_wy * pi * pi * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * mu * w_z * cos(a_wz * pi * z / L) * a_wz * a_wz * pi * pi * pow(L, -0.2e1) - p_z * sin(a_pz * pi * z / L) * a_pz * pi / L + rho_x * cos(a_rhox * pi * x / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_rhox * pi / L - rho_y * sin(a_rhoy * pi * y / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_rhoy * pi / L + rho_z * cos(a_rhoz * pi * z / L) * pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), 0.2e1) * a_rhoz * pi / L + u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * a_ux * pi / L + v_y * cos(a_vy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * a_vy * pi / L + w_x * cos(a_wx * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_wx * pi / L + w_y * cos(a_wy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_wy * pi / L - 0.2e1 * w_z * sin(a_wz * pi * z / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * a_wz * pi / L;
  return(Q_w);
}

double eval_3d_source_rho(
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
  double L,
  double R,
  double k)
{
  const double pi = acos(-1);
  double Q_rho;
  Q_rho = rho_x * cos(a_rhox * pi * x / L) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_rhox * pi / L - rho_y * sin(a_rhoy * pi * y / L) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_rhoy * pi / L + rho_z * cos(a_rhoz * pi * z / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * a_rhoz * pi / L + u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * a_ux * pi / L + v_y * cos(a_vy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * a_vy * pi / L - w_z * sin(a_wz * pi * z / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * a_wz * pi / L;
  return(Q_rho);
}
