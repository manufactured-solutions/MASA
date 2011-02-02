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
  double RHO;
  double U;
  double P;
  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  P = p_0 + p_x * cos(a_px * pi * x / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L);

  Q_e = cos(a_rhox * pi * x / L) * pow(U, 0.3e1) * a_rhox * pi * rho_x / L / 0.2e1 + cos(a_ux * pi * x / L) * P * a_ux * pi * u_x * Gamma / L / (Gamma - 0.1e1) + 0.3e1 / 0.2e1 * cos(a_ux * pi * x / L) * RHO * U * U * a_ux * pi * u_x / L - sin(a_px * pi * x / L) * U * a_px * pi * p_x * Gamma / L / (Gamma - 0.1e1);

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
  double RHO;
  double U;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L);

  Q_u = 0.2e1 * cos(a_ux * pi * x / L) * RHO * U * a_ux * pi * u_x / L + cos(a_rhox * pi * x / L) * U * U * a_rhox * pi * rho_x / L - sin(a_px * pi * x / L) * a_px * pi * p_x / L;

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
  double RHO;
  double U;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L);

  Q_rho = cos(a_ux * pi * x / L) * RHO * a_ux * pi * u_x / L + cos(a_rhox * pi * x / L) * U * a_rhox * pi * rho_x / L;

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
  double RHO;
  double U;
  double V;
  double P;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);
  P = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L);

  Q_e = -a_px * pi * p_x * Gamma * U * sin(a_px * pi * x / L) / (Gamma - 0.1e1) / L + a_py * pi * p_y * Gamma * V * cos(a_py * pi * y / L) / (Gamma - 0.1e1) / L + (U * U + V * V) * a_rhox * pi * rho_x * U * cos(a_rhox * pi * x / L) / L / 0.2e1 - (U * U + V * V) * a_rhoy * pi * rho_y * V * sin(a_rhoy * pi * y / L) / L / 0.2e1 + (0.3e1 * a_ux * u_x * cos(a_ux * pi * x / L) + a_vy * v_y * cos(a_vy * pi * y / L)) * pi * RHO * U * U / L / 0.2e1 - (a_uy * u_y * sin(a_uy * pi * y / L) + a_vx * v_x * sin(a_vx * pi * x / L)) * pi * RHO * U * V / L + (a_ux * u_x * cos(a_ux * pi * x / L) + 0.3e1 * a_vy * v_y * cos(a_vy * pi * y / L)) * pi * RHO * V * V / L / 0.2e1 + (a_ux * u_x * cos(a_ux * pi * x / L) + a_vy * v_y * cos(a_vy * pi * y / L)) * pi * Gamma * P / (Gamma - 0.1e1) / L;

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
  double RHO;
  double U;
  double V;
  
  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);

  Q_u = a_rhox * pi * rho_x * U * U * cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * U * V * sin(a_rhoy * pi * y / L) / L - a_uy * pi * u_y * RHO * V * sin(a_uy * pi * y / L) / L - a_px * pi * p_x * sin(a_px * pi * x / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * pi * x / L) + a_vy * v_y * cos(a_vy * pi * y / L)) * pi * RHO * U / L;
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
  double RHO;
  double U;
  double V;
  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);

  Q_v = a_rhox * pi * rho_x * U * V * cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * V * V * sin(a_rhoy * pi * y / L) / L - a_vx * pi * v_x * RHO * U * sin(a_vx * pi * x / L) / L + a_py * pi * p_y * cos(a_py * pi * y / L) / L + (a_ux * u_x * cos(a_ux * pi * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * pi * y / L)) * pi * RHO * V / L;
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
  double RHO;
  double U;
  double V;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);

  Q_rho = a_rhox * pi * rho_x * U * cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * V * sin(a_rhoy * pi * y / L) / L + (a_ux * u_x * cos(a_ux * pi * x / L) + a_vy * v_y * cos(a_vy * pi * y / L)) * pi * RHO / L;
  return(Q_rho);
}

//--------------------------------------------------------------------------
// 
// 
//    3D problem
// 
//
//--------------------------------------------------------------------------

double eval_3d_p_an(double x,double y,double z,double p_0,double p_x,double p_y,double p_z,double a_px,double a_py,double a_pz,double L)
{
  const double pi = acos(-1);
  double p_an = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L);
  return p_an;
}
  
double eval_3d_u_an(double x,double y,double z,double u_0,double u_x,double u_y,double u_z,double a_ux,double a_uy,double a_uz,double L)
{
  const double pi = acos(-1);
  double u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L);  
  return u_an;
} 
 
double eval_3d_v_an(double x,double y,double z,double v_0,double v_x,double v_y,double v_z,double a_vx,double a_vy,double a_vz,double L)
{
  const double pi = acos(-1);
  double v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L);
  return v_an;
}

double eval_3d_w_an(double x,double y,double z,double w_0,double w_x,double w_y,double w_z,double a_wx,double a_wy,double a_wz,double L)
{
  const double pi = acos(-1);
  double w_an = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L);
  return w_an;
}

double eval_3d_rho_an(double x,double y,double z,double rho_0,double rho_x,double rho_y,double rho_z,double a_rhox,double a_rhoy,double a_rhoz,double L)
{ 
  const double pi = acos(-1);
  double rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L);
  return rho_an;
}


double eval_3d_e_source(  
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
  double RHO;
  double P;
  double U;
  double V;
  double W;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L);
  W = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L);
  P = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L);

  Q_e = -a_px * pi * p_x * Gamma * U * sin(a_px * pi * x / L) / (Gamma - 0.1e1) / L + a_py * pi * p_y * Gamma * V * cos(a_py * pi * y / L) / (Gamma - 0.1e1) / L - a_pz * pi * p_z * Gamma * W * sin(a_pz * pi * z / L) / (Gamma - 0.1e1) / L + (U * U + V * V + W * W) * a_rhox * pi * rho_x * U * cos(a_rhox * pi * x / L) / L / 0.2e1 - (U * U + V * V + W * W) * a_rhoy * pi * rho_y * V * sin(a_rhoy * pi * y / L) / L / 0.2e1 + (U * U + V * V + W * W) * a_rhoz * pi * rho_z * W * cos(a_rhoz * pi * z / L) / L / 0.2e1 - (-0.3e1 * a_ux * u_x * cos(a_ux * pi * x / L) - a_vy * v_y * cos(a_vy * pi * y / L) + a_wz * w_z * sin(a_wz * pi * z / L)) * pi * RHO * U * U / L / 0.2e1 - (a_uy * u_y * sin(a_uy * pi * y / L) + a_vx * v_x * sin(a_vx * pi * x / L)) * pi * RHO * U * V / L - (a_uz * u_z * sin(a_uz * pi * z / L) - a_wx * w_x * cos(a_wx * pi * x / L)) * pi * RHO * U * W / L - (-a_ux * u_x * cos(a_ux * pi * x / L) - 0.3e1 * a_vy * v_y * cos(a_vy * pi * y / L) + a_wz * w_z * sin(a_wz * pi * z / L)) * pi * RHO * V * V / L / 0.2e1 + (a_vz * v_z * cos(a_vz * pi * z / L) + a_wy * w_y * cos(a_wy * pi * y / L)) * pi * RHO * V * W / L - (-a_ux * u_x * cos(a_ux * pi * x / L) - a_vy * v_y * cos(a_vy * pi * y / L) + 0.3e1 * a_wz * w_z * sin(a_wz * pi * z / L)) * pi * RHO * W * W / L / 0.2e1 - (-a_ux * u_x * cos(a_ux * pi * x / L) - a_vy * v_y * cos(a_vy * pi * y / L) + a_wz * w_z * sin(a_wz * pi * z / L)) * pi * Gamma * P / (Gamma - 0.1e1) / L;

  return(Q_e);
}

double eval_3d_u_source(
  double x,  double y,  double z,
  double u_0,  double u_x,  double u_y,  double u_z,
  double v_0,  double v_x,  double v_y,  double v_z,
  double w_0,  double w_x,  double w_y,  double w_z,
  double rho_0,  double rho_x,  double rho_y,  double rho_z,
  double p_0,  double p_x,  double p_y,  double p_z,
  double a_px,  double a_py,  double a_pz,
  double a_rhox,  double a_rhoy,  double a_rhoz,
  double a_ux,  double a_uy,  double a_uz,
  double a_vx,  double a_vy,  double a_vz,
  double a_wx,  double a_wy,  double a_wz,
  double L)
{
  const double pi = acos(-1);
  double Q_u;
  double RHO;
  double U;
  double V;
  double W;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L);
  W = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L);

  Q_u = a_rhox * pi * rho_x * U * U * cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * U * V * sin(a_rhoy * pi * y / L) / L + a_rhoz * pi * rho_z * U * W * cos(a_rhoz * pi * z / L) / L - a_uy * pi * u_y * RHO * V * sin(a_uy * pi * y / L) / L - a_uz * pi * u_z * RHO * W * sin(a_uz * pi * z / L) / L - a_px * pi * p_x * sin(a_px * pi * x / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * pi * x / L) + a_vy * v_y * cos(a_vy * pi * y / L) - a_wz * w_z * sin(a_wz * pi * z / L)) * pi * RHO * U / L;

  return(Q_u);
}

double eval_3d_v_source(
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
  double RHO;
  double U;
  double V;
  double W;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L);
  W = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L);

  Q_v = a_rhox * pi * rho_x * U * V * cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * V * V * sin(a_rhoy * pi * y / L) / L + a_rhoz * pi * rho_z * V * W * cos(a_rhoz * pi * z / L) / L - a_vx * pi * v_x * RHO * U * sin(a_vx * pi * x / L) / L + a_vz * pi * v_z * RHO * W * cos(a_vz * pi * z / L) / L + a_py * pi * p_y * cos(a_py * pi * y / L) / L + (a_ux * u_x * cos(a_ux * pi * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * pi * y / L) - a_wz * w_z * sin(a_wz * pi * z / L)) * pi * RHO * V / L;
  return(Q_v);
}

double eval_3d_w_source(
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
  double RHO;
  double U;
  double V;
  double W;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L);
  W = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L);

  Q_w = a_rhox * pi * rho_x * U * W * cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * V * W * sin(a_rhoy * pi * y / L) / L + a_rhoz * pi * rho_z * W * W * cos(a_rhoz * pi * z / L) / L + a_wx * pi * w_x * RHO * U * cos(a_wx * pi * x / L) / L + a_wy * pi * w_y * RHO * V * cos(a_wy * pi * y / L) / L - a_pz * pi * p_z * sin(a_pz * pi * z / L) / L + (a_ux * u_x * cos(a_ux * pi * x / L) + a_vy * v_y * cos(a_vy * pi * y / L) - 0.2e1 * a_wz * w_z * sin(a_wz * pi * z / L)) * pi * RHO * W / L;
  return(Q_w);
}

double eval_3d_rho_source(
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
  double RHO;
  double U;
  double V;
  double W;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L);
  W = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L);

  Q_rho = a_rhox * pi * rho_x * U * cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * V * sin(a_rhoy * pi * y / L) / L + a_rhoz * pi * rho_z * W * cos(a_rhoz * pi * z / L) / L + (a_ux * u_x * cos(a_ux * pi * x / L) + a_vy * v_y * cos(a_vy * pi * y / L) - a_wz * w_z * sin(a_wz * pi * z / L)) * pi * RHO / L;

  return(Q_rho);
}
