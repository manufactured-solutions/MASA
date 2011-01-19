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
// $Author: nick $
// $Id: euler-source.c 12693 2010-08-26 03:35:34Z nick $
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <math.h>

double eval_1d_p_an(double x,double p_0,double p_x,double a_px,double L)
{
  const double pi = acos(-1);
  double p_an = p_0 + p_x * cos(a_px * pi * x / L);
  return p_an;
}
  
double eval_1d_u_an(double x,double u_0,double u_x,double a_ux,double L)
{
  const double pi = acos(-1);
 
  double u_an = u_0 + u_x * sin(a_ux * pi * x / L);
  return u_an;
} 
 
double eval_1d_rho_an(double x,double rho_0,double rho_x,double a_rhox,double L)
{ 
  const double pi = acos(-1);

  double rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  return rho_an;
}

double eval_1d_e_source(
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
  const double pi = acos(-1);

  double Q_e;
  Q_e = cos(a_rhox * pi * x / L) * rho_x * pow(u_0 + u_x * sin(a_ux * pi * x / L), 0.3e1) * a_rhox * pi / L / 0.2e1 + cos(a_ux * pi * x / L) * (p_0 + p_x * cos(a_px * pi * x / L)) * a_ux * pi * u_x * Gamma / L / (Gamma - 0.1e1) - Gamma * p_x * sin(a_px * pi * x / L) * (u_0 + u_x * sin(a_ux * pi * x / L)) * a_px * pi / L / (Gamma - 0.1e1) + 0.3e1 / 0.2e1 * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L)) * pow(u_0 + u_x * sin(a_ux * pi * x / L), 0.2e1) * a_ux * pi * u_x / L;
  return(Q_e);
}

double eval_1d_u_source(
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
  const double pi = acos(-1);

  double Q_u;
  Q_u = -sin(a_px * pi * x / L) * a_px * pi * p_x / L + rho_x * cos(a_rhox * pi * x / L) * pow(u_0 + u_x * sin(a_ux * pi * x / L), 0.2e1) * a_rhox * pi / L + 0.2e1 * u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L)) * (u_0 + u_x * sin(a_ux * pi * x / L)) * a_ux * pi / L;
  return(Q_u);
}

double eval_1d_rho_source(  
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
  const double pi = acos(-1);
 
  double Q_rho;
  Q_rho = rho_x * cos(a_rhox * pi * x / L) * (u_0 + u_x * sin(a_ux * pi * x / L)) * a_rhox * pi / L + u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L)) * a_ux * pi / L;
  return(Q_rho);
}


//--------------------------------------------------------------------------
// 
// 
//    2D problem
// 
//
//--------------------------------------------------------------------------

double eval_2d_p_an(double x,double y,double p_0,double p_x,double p_y,double a_px,double a_py,double L)
{
  const double pi = acos(-1);
  double p_an = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L);
  return p_an;
}
  
double eval_2d_u_an(double x,double y,double u_0,double u_x,double u_y,double a_ux,double a_uy,double L)
{
  const double pi = acos(-1);
  double u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  return u_an;
} 
 
double eval_2d_v_an(double x,double y,double v_0,double v_x,double v_y,double a_vx,double a_vy,double L)
{
  const double pi = acos(-1);
  double v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);
  return v_an;
}

double eval_2d_rho_an(double x,double y,double rho_0,double rho_x,double rho_y,double a_rhox,double a_rhoy,double L)
{ 
  const double pi = acos(-1);
  double rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  return rho_an;
}

double eval_2d_e_source(
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
			 double L)
{
  const double pi = acos(-1);
  double Q_e;
  Q_e = -Gamma * (u_x * sin(x * a_ux * pi / L) + u_y * cos(a_uy * pi * y / L) + u_0) * a_px * p_x * pi * sin(a_px * pi * x / L) / (Gamma - 0.1e1) / L + Gamma * (v_x * cos(a_vx * pi * x / L) + v_y * sin(y * a_vy * pi / L) + v_0) * a_py * p_y * pi * cos(a_py * pi * y / L) / (Gamma - 0.1e1) / L + a_rhox * pi * rho_x * cos(a_rhox * pi * x / L) * (u_x * sin(x * a_ux * pi / L) + u_y * cos(a_uy * pi * y / L) + u_0) * (pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(y * a_vy * pi / L) + v_0, 0.2e1) + pow(u_x * sin(x * a_ux * pi / L) + u_y * cos(a_uy * pi * y / L) + u_0, 0.2e1)) / L / 0.2e1 - a_rhoy * pi * rho_y * sin(a_rhoy * pi * y / L) * (v_x * cos(a_vx * pi * x / L) + v_y * sin(y * a_vy * pi / L) + v_0) * (pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(y * a_vy * pi / L) + v_0, 0.2e1) + pow(u_x * sin(x * a_ux * pi / L) + u_y * cos(a_uy * pi * y / L) + u_0, 0.2e1)) / L / 0.2e1 + (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_0) * a_ux * pi * u_x * cos(x * a_ux * pi / L) * Gamma / (Gamma - 0.1e1) / L + (pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(y * a_vy * pi / L) + v_0, 0.2e1) + 0.3e1 * pow(u_x * sin(x * a_ux * pi / L) + u_y * cos(a_uy * pi * y / L) + u_0, 0.2e1)) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * a_ux * pi * u_x * cos(x * a_ux * pi / L) / L / 0.2e1 + (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_0) * a_vy * pi * v_y * cos(y * a_vy * pi / L) * Gamma / (Gamma - 0.1e1) / L + (0.3e1 * pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(y * a_vy * pi / L) + v_0, 0.2e1) + pow(u_x * sin(x * a_ux * pi / L) + u_y * cos(a_uy * pi * y / L) + u_0, 0.2e1)) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * a_vy * pi * v_y * cos(y * a_vy * pi / L) / L / 0.2e1 - (v_x * cos(a_vx * pi * x / L) + v_y * sin(y * a_vy * pi / L) + v_0) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * (u_x * sin(x * a_ux * pi / L) + u_y * cos(a_uy * pi * y / L) + u_0) * pi * a_uy * u_y * sin(a_uy * pi * y / L) / L - (v_x * cos(a_vx * pi * x / L) + v_y * sin(y * a_vy * pi / L) + v_0) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * (u_x * sin(x * a_ux * pi / L) + u_y * cos(a_uy * pi * y / L) + u_0) * pi * a_vx * v_x * sin(a_vx * pi * x / L) / L;
  return(Q_e);
}

double eval_2d_u_source(
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
			 double L
			 )
{
  const double pi = acos(-1);
  double Q_u;
  Q_u = -p_x * sin(a_px * pi * x / L) * a_px * pi / L + rho_x * cos(a_rhox * pi * x / L) * pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L), 0.2e1) * a_rhox * pi / L - rho_y * sin(a_rhoy * pi * y / L) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_rhoy * pi / L + 0.2e1 * u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_ux * pi / L - u_y * sin(a_uy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * a_uy * pi / L + v_y * cos(a_vy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_vy * pi / L;
  return(Q_u);
}

double eval_2d_v_source(
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
			 double L)
{
  const double pi = acos(-1);
  double Q_v;
  Q_v = p_y * cos(a_py * pi * y / L) * a_py * pi / L + rho_x * cos(a_rhox * pi * x / L) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_rhox * pi / L - sin(a_rhoy * pi * y / L) * rho_y * pow(v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L), 0.2e1) * a_rhoy * pi / L + cos(a_ux * pi * x / L) * u_x * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * a_ux * pi / L - sin(a_vx * pi * x / L) * v_x * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_vx * pi / L + 0.2e1 * cos(a_vy * pi * y / L) * v_y * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * a_vy * pi / L;
  return(Q_v);
}

double eval_2d_rho_source(
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
			   double L)
{
  const double pi = acos(-1);
  double Q_rho;
  Q_rho = rho_x * cos(a_rhox * pi * x / L) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_rhox * pi / L - rho_y * sin(a_rhoy * pi * y / L) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * a_rhoy * pi / L + u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * a_ux * pi / L + v_y * cos(a_vy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * a_vy * pi / L;
  return(Q_rho);
}

//--------------------------------------------------------------------------
// 
// 
//    3D problem
// 
//
//--------------------------------------------------------------------------

double eval_3d_p_an_(double x,double y,double z,double p_0,double p_x,double p_y,double p_z,double a_px,double a_py,double a_pz,double L)
{
  const double pi = acos(-1);
  double p_an = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L);
  return p_an;
}
  
double eval_3d_u_an_(double x,double y,double z,double u_0,double u_x,double u_y,double u_z,double a_ux,double a_uy,double a_uz,double L)
{
  const double pi = acos(-1);
  double u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L);  
  return u_an;
} 
 
double eval_3d_v_an_(double x,double y,double z,double v_0,double v_x,double v_y,double v_z,double a_vx,double a_vy,double a_vz,double L)
{
  const double pi = acos(-1);
  double v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L);
  return v_an;
}

double eval_3d_w_an_(double x,double y,double z,double w_0,double w_x,double w_y,double w_z,double a_wx,double a_wy,double a_wz,double L)
{
  const double pi = acos(-1);
  double w_an = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L);
  return w_an;
}

double eval_3d_rho_an_(double x,double y,double z,double rho_0,double rho_x,double rho_y,double rho_z,double a_rhox,double a_rhoy,double a_rhoz,double L)
{ 
  const double pi = acos(-1);
  double rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L);
  return rho_an;
}


double eval_3d_e_source_(  
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
  const double pi = acos(-1);
  double Q_e;
  Q_e = -Gamma * p_x * sin(a_px * pi * x / L) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_px * pi / L / (Gamma - 0.1e1) + Gamma * p_y * cos(a_py * pi * y / L) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_py * pi / L / (Gamma - 0.1e1) - Gamma * p_z * sin(a_pz * pi * z / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * a_pz * pi / L / (Gamma - 0.1e1) + (pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), 0.2e1) + pow(v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), 0.2e1)) * a_rhox * pi * rho_x * cos(a_rhox * pi * x / L) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) / L / 0.2e1 - (pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), 0.2e1) + pow(v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), 0.2e1)) * a_rhoy * pi * rho_y * sin(a_rhoy * pi * y / L) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) / L / 0.2e1 + (pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), 0.2e1) + pow(v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), 0.2e1)) * a_rhoz * pi * rho_z * cos(a_rhoz * pi * z / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) / L / 0.2e1 + a_ux * pi * u_x * cos(a_ux * pi * x / L) * ((0.3e1 * rho_x * sin(a_rhox * pi * x / L) + 0.3e1 * rho_y * cos(a_rhoy * pi * y / L) + 0.3e1 * rho_z * sin(a_rhoz * pi * z / L) + 0.3e1 * rho_0) * pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), 0.2e1) + (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_0) * pow(v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L), 0.2e1) + (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_0) * pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), 0.2e1)) / L / 0.2e1 + Gamma * a_ux * pi * u_x * cos(a_ux * pi * x / L) * (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_0) / L / (Gamma - 0.1e1) - (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_0) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * pi * a_uy * u_y * sin(a_uy * pi * y / L) / L - (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_0) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * pi * a_uz * u_z * sin(a_uz * pi * z / L) / L - (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_0) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * pi * a_vx * v_x * sin(a_vx * pi * x / L) / L + a_vy * pi * v_y * cos(a_vy * pi * y / L) * ((rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_0) * pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), 0.2e1) + (0.3e1 * rho_x * sin(a_rhox * pi * x / L) + 0.3e1 * rho_y * cos(a_rhoy * pi * y / L) + 0.3e1 * rho_z * sin(a_rhoz * pi * z / L) + 0.3e1 * rho_0) * pow(v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L), 0.2e1) + (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_0) * pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), 0.2e1)) / L / 0.2e1 + Gamma * a_vy * pi * v_y * cos(a_vy * pi * y / L) * (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_0) / L / (Gamma - 0.1e1) + (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_0) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * pi * a_vz * v_z * cos(a_vz * pi * z / L) / L + (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_0) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * pi * a_wx * w_x * cos(a_wx * pi * x / L) / L + (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_0) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * pi * a_wy * w_y * cos(a_wy * pi * y / L) / L - a_wz * pi * w_z * sin(a_wz * pi * z / L) * ((rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_0) * pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), 0.2e1) + (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_0) * pow(v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L), 0.2e1) + (0.3e1 * rho_x * sin(a_rhox * pi * x / L) + 0.3e1 * rho_y * cos(a_rhoy * pi * y / L) + 0.3e1 * rho_z * sin(a_rhoz * pi * z / L) + 0.3e1 * rho_0) * pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), 0.2e1)) / L / 0.2e1 - Gamma * a_wz * pi * w_z * sin(a_wz * pi * z / L) * (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_0) / L / (Gamma - 0.1e1);
  return(Q_e);
}

double eval_3d_u_source_(
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
  const double pi = acos(-1);
  double Q_u;
  Q_u = -p_x * sin(a_px * pi * x / L) * a_px * pi / L + rho_x * cos(a_rhox * pi * x / L) * pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), 0.2e1) * a_rhox * pi / L - rho_y * sin(a_rhoy * pi * y / L) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_rhoy * pi / L + rho_z * cos(a_rhoz * pi * z / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_rhoz * pi / L + 0.2e1 * u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_ux * pi / L - u_y * sin(a_uy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_uy * pi / L - u_z * sin(a_uz * pi * z / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * a_uz * pi / L + v_y * cos(a_vy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_vy * pi / L - w_z * sin(a_wz * pi * z / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_wz * pi / L;
  return(Q_u);
}

double eval_3d_v_source_(
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
  const double pi = acos(-1);
  double Q_v;
  Q_v = p_y * cos(a_py * pi * y / L) * a_py * pi / L + rho_x * cos(a_rhox * pi * x / L) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_rhox * pi / L - rho_y * sin(a_rhoy * pi * y / L) * pow(v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L), 0.2e1) * a_rhoy * pi / L + rho_z * cos(a_rhoz * pi * z / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_rhoz * pi / L + u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_ux * pi / L - v_x * sin(a_vx * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_vx * pi / L + 0.2e1 * v_y * cos(a_vy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_vy * pi / L + v_z * cos(a_vz * pi * z / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * a_vz * pi / L - w_z * sin(a_wz * pi * z / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_wz * pi / L;
  return(Q_v);
}

double eval_3d_w_source_(
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
  const double pi = acos(-1);
  double Q_w;
  Q_w = -p_z * sin(a_pz * pi * z / L) * a_pz * pi / L + rho_x * cos(a_rhox * pi * x / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_rhox * pi / L - rho_y * sin(a_rhoy * pi * y / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_rhoy * pi / L + rho_z * cos(a_rhoz * pi * z / L) * pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), 0.2e1) * a_rhoz * pi / L + u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * a_ux * pi / L + v_y * cos(a_vy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * a_vy * pi / L + w_x * cos(a_wx * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_wx * pi / L + w_y * cos(a_wy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_wy * pi / L - 0.2e1 * w_z * sin(a_wz * pi * z / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * a_wz * pi / L;
  return(Q_w);
}

double eval_3d_rho_source_(
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

  const double pi = acos(-1);
  double Q_rho;
  Q_rho = rho_x * cos(a_rhox * pi * x / L) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_rhox * pi / L - rho_y * sin(a_rhoy * pi * y / L) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_rhoy * pi / L + rho_z * cos(a_rhoz * pi * z / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * a_rhoz * pi / L + u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * a_ux * pi / L + v_y * cos(a_vy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * a_vy * pi / L - w_z * sin(a_wz * pi * z / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * a_wz * pi / L;
  return(Q_rho);
}
