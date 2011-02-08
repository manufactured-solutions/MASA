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

#ifdef MASA_STRICT_REGRESSION
const double threshold = 1.0e-15; // should be small enough to catch any obvious problems
#else
// forward error propagation: worst case scenario FP error is infinity, wouldn't expect FP error of less than 1000*epsilon. 
// My best guess of the error level would have been log2(1000)*epsilon, or roughly 10*epsilon. Seeing 8*epsilon. 
// Thus, weakening for distribution 
const double threshold = 2.0e-15; // weakening because of O(1000) FP operations failing on some macbooks
#endif

double anQ_p (double x,double y,double z,double p_0,double p_x,double p_y,double p_z,double a_px,double a_py,double a_pz,double L)
{
  double pi = acos(-1);
  double exact_p = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L);
  return exact_p;
}
 
double anQ_u (double x,double y,double z,double u_0,double u_x,double u_y,double u_z,double a_ux,double a_uy,double a_uz,double L)
{
  double pi = acos(-1); 
  double exact_u = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L);  
  return exact_u;
} 
 
double anQ_v (double x,double y,double z,double v_0,double v_x,double v_y,double v_z,double a_vx,double a_vy,double a_vz,double L)
{
  double pi = acos(-1);
  double exact_v = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L);
  return exact_v;
}

double anQ_w (double x,double y,double z,double w_0,double w_x,double w_y,double w_z,double a_wx,double a_wy,double a_wz,double L)
{
  double pi = acos(-1);
  double exact_w = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L);
  return exact_w;
}

double anQ_rho (double x,double y,double z,double rho_0,double rho_x,double rho_y,double rho_z,double a_rhox,double a_rhoy,double a_rhoz,double L)
{ 
  double pi = acos(-1);
  double exact_rho = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L);
  return exact_rho;
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
  double pi = acos(-1);
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
  double pi = acos(-1);
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
  double pi = acos(-1);
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
  double pi = acos(-1);
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
  double pi = acos(-1);
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

  double exact_u,exact_u2,exact_u3;
  double exact_v,exact_v2,exact_v3;
  double exact_w,exact_w2,exact_w3;
  double exact_p,exact_p2,exact_p3;
  double exact_rho,exact_rho2,exact_rho3;

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
	  ufield = cmasa_eval_3d_source_rho_u (x,y,z);
	  vfield = cmasa_eval_3d_source_rho_v (x,y,z);
	  wfield = cmasa_eval_3d_source_rho_w (x,y,z);
	  efield = cmasa_eval_3d_source_e     (x,y,z);
	  rho    = cmasa_eval_3d_source_rho   (x,y,z);
	  
	  //evaluate analytical terms
	  exact_u   = cmasa_eval_3d_exact_u      (x,y,z);
	  exact_v   = cmasa_eval_3d_exact_v      (x,y,z);
	  exact_w   = cmasa_eval_3d_exact_w      (x,y,z);
	  exact_p   = cmasa_eval_3d_exact_p      (x,y,z);
	  exact_rho = cmasa_eval_3d_exact_rho     (x,y,z);

	  // check against maple output
	  ufield2   = SourceQ_u  (x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,L);
	  vfield2   = SourceQ_v  (x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,L);
	  wfield2   = SourceQ_w  (x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,L);
	  rho2      = SourceQ_rho(x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,mu,L);
	  efield2   = SourceQ_e  (x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,mu,Gamma,L);
  
	  exact_u2     = anQ_u   (x,y,z,u_0,u_x,u_y,u_z,a_ux,a_uy,a_uz,L);
	  exact_v2     = anQ_v   (x,y,z,v_0,v_x,v_y,v_z,a_vx,a_vy,a_vz,L);
	  exact_w2     = anQ_w   (x,y,z,w_0,w_x,w_y,w_z,a_wx,a_wy,a_wz,L);
	  exact_rho2   = anQ_rho (x,y,z,rho_0,rho_x,rho_y,rho_z,a_rhox,a_rhoy,a_rhoz,L);
	  exact_p2     = anQ_p   (x,y,z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,L);

	  // test the result is roughly zero
	  // choose between abs and rel error
#ifdef MASA_STRICT_REGRESSION

	  ufield3 = fabs(ufield-ufield2);
	  vfield3 = fabs(vfield-vfield2);
	  wfield3 = fabs(wfield-wfield2);
	  efield3 = fabs(efield-efield2);
	  rho3    = fabs(rho-rho2);

	  exact_u3   = fabs(exact_u-exact_u2);
	  exact_v3   = fabs(exact_v-exact_v2);
	  exact_w3   = fabs(exact_w-exact_w2);
	  exact_rho3 = fabs(exact_rho-exact_rho2);
	  exact_p3   = fabs(exact_p-exact_p2);

#else

	  ufield3 = fabs(ufield-ufield2)/fabs(ufield2);
	  vfield3 = fabs(vfield-vfield2)/fabs(vfield2);
	  wfield3 = fabs(wfield-wfield2)/fabs(wfield2);
	  efield3 = fabs(efield-efield2)/fabs(efield2);
	  rho3    = fabs(rho-rho2)/fabs(rho2);

	  exact_u3   = fabs(exact_u-exact_u2)/fabs(exact_u2);
	  exact_v3   = fabs(exact_v-exact_v2)/fabs(exact_v2);
	  exact_w3   = fabs(exact_w-exact_w2)/fabs(exact_w2);
	  exact_rho3 = fabs(exact_rho-exact_rho2)/fabs(exact_rho2);
	  exact_p3   = fabs(exact_p-exact_p2)/fabs(exact_p2);

#endif

	  if(ufield3 > threshold)
	    {
	      printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-3d\n");
	      printf("U Field Source Term\n");
	      exit(1);
	    }
	  
	  if(exact_u3 > threshold)
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

	  if(exact_v3 > threshold)
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

	  if(exact_w3 > threshold)
	    {
	      printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-3d\n");
	      printf("W Field Analytical Term\n");
	      printf("Threshold Exceeded: %g\n\n",exact_w);
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

	  if(exact_p3 > threshold)
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

	  if(exact_rho3 > threshold)
	    {
	      printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-3d\n");
	      printf("RHO Field Analytical Term\n");
	      exit(1);
	    }
	  
	}// done iterating

  // tests passed
  return 0;
}
