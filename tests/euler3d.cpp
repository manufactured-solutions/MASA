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
// euler3d.cpp: program that tests euler3d from masa against known source
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <limits>

using namespace MASA;
using namespace std;

typedef double Scalar;

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

template<typename Scalar>
Scalar anQ_p (Scalar x,Scalar y,Scalar z,Scalar p_0,Scalar p_x,Scalar p_y,Scalar p_z,Scalar a_px,Scalar a_py,Scalar a_pz,Scalar L)
{
  Scalar pi = acos(-1);
  Scalar exact_p = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L);
  return exact_p;
}
 
template<typename Scalar>
Scalar anQ_u (Scalar x,Scalar y,Scalar z,Scalar u_0,Scalar u_x,Scalar u_y,Scalar u_z,Scalar a_ux,Scalar a_uy,Scalar a_uz,Scalar L)
{
  Scalar pi = acos(-1);
  Scalar exact_u = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L);  
  return exact_u;
} 

 template<typename Scalar>
Scalar anQ_v (Scalar x,Scalar y,Scalar z,Scalar v_0,Scalar v_x,Scalar v_y,Scalar v_z,Scalar a_vx,Scalar a_vy,Scalar a_vz,Scalar L)
{
  Scalar pi = acos(-1);
  Scalar exact_v = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L);
  return exact_v;
}

template<typename Scalar>
Scalar anQ_w (Scalar x,Scalar y,Scalar z,Scalar w_0,Scalar w_x,Scalar w_y,Scalar w_z,Scalar a_wx,Scalar a_wy,Scalar a_wz,Scalar L)
{
  Scalar pi = acos(-1);
  Scalar exact_w = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L);
  return exact_w;
}

template<typename Scalar>
Scalar anQ_rho (Scalar x,Scalar y,Scalar z,Scalar rho_0,Scalar rho_x,Scalar rho_y,Scalar rho_z,Scalar a_rhox,Scalar a_rhoy,Scalar a_rhoz,Scalar L)
{ 
  Scalar pi = acos(-1);
  Scalar exact_rho = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L);
  return exact_rho;
}

template<typename Scalar>
Scalar SourceQ_e (
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar u_z,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar v_z,
  Scalar w_0,
  Scalar w_x,
  Scalar w_y,
  Scalar w_z,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar rho_z,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar p_z,
  Scalar a_px,
  Scalar a_py,
  Scalar a_pz,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_rhoz,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_uz,
  Scalar a_vx,
  Scalar a_vy,
  Scalar a_vz,
  Scalar a_wx,
  Scalar a_wy,
  Scalar a_wz,
  Scalar mu,
  Scalar Gamma,
  Scalar L)
{
  Scalar pi = acos(-1);
  Scalar Q_e;
  Scalar RHO;
  Scalar P;
  Scalar U;
  Scalar V;
  Scalar W;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L);
  W = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L);
  P = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L);

  Q_e = -a_px * pi * p_x * Gamma * U * sin(a_px * pi * x / L) / (Gamma - 0.1e1) / L + a_py * pi * p_y * Gamma * V * cos(a_py * pi * y / L) / (Gamma - 0.1e1) / L - a_pz * pi * p_z * Gamma * W * sin(a_pz * pi * z / L) / (Gamma - 0.1e1) / L + (U * U + V * V + W * W) * a_rhox * pi * rho_x * U * cos(a_rhox * pi * x / L) / L / 0.2e1 - (U * U + V * V + W * W) * a_rhoy * pi * rho_y * V * sin(a_rhoy * pi * y / L) / L / 0.2e1 + (U * U + V * V + W * W) * a_rhoz * pi * rho_z * W * cos(a_rhoz * pi * z / L) / L / 0.2e1 - (-0.3e1 * a_ux * u_x * cos(a_ux * pi * x / L) - a_vy * v_y * cos(a_vy * pi * y / L) + a_wz * w_z * sin(a_wz * pi * z / L)) * pi * RHO * U * U / L / 0.2e1 - (a_uy * u_y * sin(a_uy * pi * y / L) + a_vx * v_x * sin(a_vx * pi * x / L)) * pi * RHO * U * V / L - (a_uz * u_z * sin(a_uz * pi * z / L) - a_wx * w_x * cos(a_wx * pi * x / L)) * pi * RHO * U * W / L - (-a_ux * u_x * cos(a_ux * pi * x / L) - 0.3e1 * a_vy * v_y * cos(a_vy * pi * y / L) + a_wz * w_z * sin(a_wz * pi * z / L)) * pi * RHO * V * V / L / 0.2e1 + (a_vz * v_z * cos(a_vz * pi * z / L) + a_wy * w_y * cos(a_wy * pi * y / L)) * pi * RHO * V * W / L - (-a_ux * u_x * cos(a_ux * pi * x / L) - a_vy * v_y * cos(a_vy * pi * y / L) + 0.3e1 * a_wz * w_z * sin(a_wz * pi * z / L)) * pi * RHO * W * W / L / 0.2e1 - (-a_ux * u_x * cos(a_ux * pi * x / L) - a_vy * v_y * cos(a_vy * pi * y / L) + a_wz * w_z * sin(a_wz * pi * z / L)) * pi * Gamma * P / (Gamma - 0.1e1) / L;

  return(Q_e);
}


template<typename Scalar>
Scalar SourceQ_u ( // should be 38
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar u_z,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar v_z,
  Scalar w_0,
  Scalar w_x,
  Scalar w_y,
  Scalar w_z,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar rho_z,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar p_z,
  Scalar a_px,
  Scalar a_py,
  Scalar a_pz,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_rhoz,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_uz,
  Scalar a_vx,
  Scalar a_vy,
  Scalar a_vz,
  Scalar a_wx,
  Scalar a_wy,
  Scalar a_wz,
  Scalar L)
{
  Scalar pi = acos(-1);
  Scalar Q_u;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar W;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L);
  W = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L);

  Q_u = a_rhox * pi * rho_x * U * U * cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * U * V * sin(a_rhoy * pi * y / L) / L + a_rhoz * pi * rho_z * U * W * cos(a_rhoz * pi * z / L) / L - a_uy * pi * u_y * RHO * V * sin(a_uy * pi * y / L) / L - a_uz * pi * u_z * RHO * W * sin(a_uz * pi * z / L) / L - a_px * pi * p_x * sin(a_px * pi * x / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * pi * x / L) + a_vy * v_y * cos(a_vy * pi * y / L) - a_wz * w_z * sin(a_wz * pi * z / L)) * pi * RHO * U / L;

  return(Q_u);
}

template<typename Scalar>
Scalar SourceQ_v ( // 38
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar u_z,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar v_z,
  Scalar w_0,
  Scalar w_x,
  Scalar w_y,
  Scalar w_z,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar rho_z,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar p_z,
  Scalar a_px,
  Scalar a_py,
  Scalar a_pz,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_rhoz,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_uz,
  Scalar a_vx,
  Scalar a_vy,
  Scalar a_vz,
  Scalar a_wx,
  Scalar a_wy,
  Scalar a_wz,
  Scalar L)
{

  Scalar pi = acos(-1);
  Scalar Q_v;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar W;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L);
  W = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L);

  Q_v = a_rhox * pi * rho_x * U * V * cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * V * V * sin(a_rhoy * pi * y / L) / L + a_rhoz * pi * rho_z * V * W * cos(a_rhoz * pi * z / L) / L - a_vx * pi * v_x * RHO * U * sin(a_vx * pi * x / L) / L + a_vz * pi * v_z * RHO * W * cos(a_vz * pi * z / L) / L + a_py * pi * p_y * cos(a_py * pi * y / L) / L + (a_ux * u_x * cos(a_ux * pi * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * pi * y / L) - a_wz * w_z * sin(a_wz * pi * z / L)) * pi * RHO * V / L;

  return(Q_v);
}

Scalar SourceQ_w ( // 38
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar u_z,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar v_z,
  Scalar w_0,
  Scalar w_x,
  Scalar w_y,
  Scalar w_z,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar rho_z,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar p_z,
  Scalar a_px,
  Scalar a_py,
  Scalar a_pz,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_rhoz,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_uz,
  Scalar a_vx,
  Scalar a_vy,
  Scalar a_vz,
  Scalar a_wx,
  Scalar a_wy,
  Scalar a_wz,
  Scalar L)
{
  Scalar pi = acos(-1);
  Scalar Q_w;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar W;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L);
  W = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L);

  Q_w = a_rhox * pi * rho_x * U * W * cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * V * W * sin(a_rhoy * pi * y / L) / L + a_rhoz * pi * rho_z * W * W * cos(a_rhoz * pi * z / L) / L + a_wx * pi * w_x * RHO * U * cos(a_wx * pi * x / L) / L + a_wy * pi * w_y * RHO * V * cos(a_wy * pi * y / L) / L - a_pz * pi * p_z * sin(a_pz * pi * z / L) / L + (a_ux * u_x * cos(a_ux * pi * x / L) + a_vy * v_y * cos(a_vy * pi * y / L) - 0.2e1 * a_wz * w_z * sin(a_wz * pi * z / L)) * pi * RHO * W / L;

  return(Q_w);
}

template<typename Scalar>
Scalar SourceQ_rho( // 39
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar u_z,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar v_z,
  Scalar w_0,
  Scalar w_x,
  Scalar w_y,
  Scalar w_z,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar rho_z,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar p_z,
  Scalar a_px,
  Scalar a_py,
  Scalar a_pz,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_rhoz,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_uz,
  Scalar a_vx,
  Scalar a_vy,
  Scalar a_vz,
  Scalar a_wx,
  Scalar a_wy,
  Scalar a_wz,
  Scalar mu,
  Scalar L)
{
  Scalar pi = acos(-1);
  Scalar Q_rho;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar W;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L);
  W = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L);

  Q_rho = a_rhox * pi * rho_x * U * cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * V * sin(a_rhoy * pi * y / L) / L + a_rhoz * pi * rho_z * W * cos(a_rhoz * pi * z / L) / L + (a_ux * u_x * cos(a_ux * pi * x / L) + a_vy * v_y * cos(a_vy * pi * y / L) - a_wz * w_z * sin(a_wz * pi * z / L)) * pi * RHO / L;

  return(Q_rho);
}

template<typename Scalar>
int run_regression()
{

  const Scalar threshold = 5 * numeric_limits<Scalar>::epsilon();

  //variables
  Scalar u_0;
  Scalar u_x;
  Scalar u_y;
  Scalar u_z;
  Scalar v_0;
  Scalar v_x;
  Scalar v_y;
  Scalar v_z;
  Scalar w_0;
  Scalar w_x;
  Scalar w_y;
  Scalar w_z;
  Scalar rho_0;
  Scalar rho_x;
  Scalar rho_y;
  Scalar rho_z;
  Scalar p_0;
  Scalar p_x;
  Scalar p_y;
  Scalar p_z;
  Scalar a_px;
  Scalar a_py;
  Scalar a_pz;
  Scalar a_rhox;
  Scalar a_rhoy;
  Scalar a_rhoz;
  Scalar a_ux;
  Scalar a_uy;
  Scalar a_uz;
  Scalar a_vx;
  Scalar a_vy;
  Scalar a_vz;
  Scalar a_wx;
  Scalar a_wy;
  Scalar a_wz;
  Scalar mu;
  Scalar Gamma;
  Scalar L;    

  // parameters
  Scalar x;
  Scalar y;
  Scalar z;

  // solutions
  Scalar ufield,ufield2,ufield3;
  Scalar vfield,vfield2,vfield3;
  Scalar wfield,wfield2,wfield3;
  Scalar efield,efield2,efield3;
  Scalar rho,rho2,rho3;
  Scalar gradx,grady,gradz;

  Scalar exact_u,exact_u2,exact_u3;
  Scalar exact_v,exact_v2,exact_v3;
  Scalar exact_w,exact_w2,exact_w3;
  Scalar exact_p,exact_p2,exact_p3;
  Scalar exact_rho,exact_rho2,exact_rho3;

  // initalize
  int nx = 33;             // number of points
  int ny = 12;  
  int nz = 45;  
  int lx=1;                // length
  int ly=2; 
  int lz=3;   
  Scalar dx=Scalar(lx)/Scalar(nx); // spacing
  Scalar dy=Scalar(ly)/Scalar(ny);
  Scalar dz=Scalar(lz)/Scalar(nz);

  masa_init<Scalar>("euler-test","euler_3d");

  // set params
  masa_init_param<Scalar>();

  // get vars for comparison
  u_0 = masa_get_param<Scalar>("u_0");
  u_x = masa_get_param<Scalar>("u_x");
  u_y = masa_get_param<Scalar>("u_y");
  u_z = masa_get_param<Scalar>("u_z");

  v_0 = masa_get_param<Scalar>("v_0");
  v_x = masa_get_param<Scalar>("v_x");
  v_y = masa_get_param<Scalar>("v_y");
  v_z = masa_get_param<Scalar>("v_z");

  w_0 = masa_get_param<Scalar>("w_0");
  w_x = masa_get_param<Scalar>("w_x");
  w_y = masa_get_param<Scalar>("w_y");
  w_z = masa_get_param<Scalar>("w_z");

  rho_0 = masa_get_param<Scalar>("rho_0");
  rho_x = masa_get_param<Scalar>("rho_x");
  rho_y = masa_get_param<Scalar>("rho_y");
  rho_z = masa_get_param<Scalar>("rho_z");

  p_0 = masa_get_param<Scalar>("p_0");
  p_x = masa_get_param<Scalar>("p_x");
  p_y = masa_get_param<Scalar>("p_y");
  p_z = masa_get_param<Scalar>("p_z");

  a_px = masa_get_param<Scalar>("a_px");
  a_py = masa_get_param<Scalar>("a_py");
  a_pz = masa_get_param<Scalar>("a_pz");

  a_rhox = masa_get_param<Scalar>("a_rhox");
  a_rhoy = masa_get_param<Scalar>("a_rhoy");
  a_rhoz = masa_get_param<Scalar>("a_rhoz");

  a_ux = masa_get_param<Scalar>("a_ux");
  a_uy = masa_get_param<Scalar>("a_uy");
  a_uz = masa_get_param<Scalar>("a_uz");

  a_vx = masa_get_param<Scalar>("a_vx");
  a_vy = masa_get_param<Scalar>("a_vy");
  a_vz = masa_get_param<Scalar>("a_vz");

  a_wx = masa_get_param<Scalar>("a_wx");
  a_wy = masa_get_param<Scalar>("a_wy");
  a_wz = masa_get_param<Scalar>("a_wz");

  Gamma = masa_get_param<Scalar>("Gamma");
  mu    = masa_get_param<Scalar>("mu");
  L     = masa_get_param<Scalar>("L");

  // check all vars initialized
  int err = masa_sanity_check<Scalar>();
  if(err != 0)
    {
      cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }

  // evaluate source terms (3D)
  for(int i=0;i<nx;i++)
    for(int j=0;j<ny;j++)    
      for(int k=0;k<nz;k++)
	{
	  x=i*dx;
	  y=j*dy;
	  z=k*dz;

	  //evalulate source terms
	  ufield = masa_eval_source_rho_u  <Scalar>(x,y,z);
	  vfield = masa_eval_source_rho_v  <Scalar>(x,y,z);
	  wfield = masa_eval_source_rho_w  <Scalar>(x,y,z);
	  efield = masa_eval_source_rho_e      <Scalar>(x,y,z);
	  rho    = masa_eval_source_rho    <Scalar>(x,y,z);
	  
	  //evaluate analytical terms
	  exact_u = masa_eval_exact_u        <Scalar>(x,y,z);
	  exact_v = masa_eval_exact_v        <Scalar>(x,y,z);
	  exact_w = masa_eval_exact_w        <Scalar>(x,y,z);
	  exact_p = masa_eval_exact_p        <Scalar>(x,y,z);
	  exact_rho = masa_eval_exact_rho    <Scalar>(x,y,z);

	  // eval gradient terms
	  gradx = masa_eval_grad_u<Scalar>(x,y,z,1);
	  grady = masa_eval_grad_u<Scalar>(x,y,z,2);
	  gradz = masa_eval_grad_u<Scalar>(x,y,z,3);

	  gradx = masa_eval_grad_v<Scalar>(x,y,z,1);
	  grady = masa_eval_grad_v<Scalar>(x,y,z,2);
	  gradz = masa_eval_grad_v<Scalar>(x,y,z,3);

	  gradx = masa_eval_grad_w<Scalar>(x,y,z,1);
	  grady = masa_eval_grad_w<Scalar>(x,y,z,2);
	  gradz = masa_eval_grad_w<Scalar>(x,y,z,3);

	  gradx = masa_eval_grad_p<Scalar>(x,y,z,1);
	  grady = masa_eval_grad_p<Scalar>(x,y,z,2);
	  gradz = masa_eval_grad_p<Scalar>(x,y,z,3);

	  gradx = masa_eval_grad_rho<Scalar>(x,y,z,1);
	  grady = masa_eval_grad_rho<Scalar>(x,y,z,2);
	  gradz = masa_eval_grad_rho<Scalar>(x,y,z,3);

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

	  nancheck(ufield3);
	  nancheck(vfield3);
	  nancheck(wfield3);
	  nancheck(efield3);
	  nancheck(rho3);
	  
	  nancheck(exact_u3);
	  nancheck(exact_v3);
	  nancheck(exact_w3);
	  nancheck(exact_rho3);
	  nancheck(exact_p3);

	  if(ufield3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "U Field Source Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << ufield3 << endl;
	      cout << "Source term is:                   " << ufield2 << endl;
	      cout << "MASA term is:                     " << ufield << endl;
	      exit(1);
	    }
	  
	  if(exact_u3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "U Field Analytical Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << exact_u << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(vfield3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "V Field Source Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << vfield3 << endl;
	      cout << "Source term is:                   " << vfield2 << endl;
	      cout << "MASA term is:                     " << vfield << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(exact_v3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "V Field Analytical Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << exact_v << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(wfield3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "W Field Source Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << wfield << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(exact_w3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "W Field Analytical Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << exact_w << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(efield3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "Energy Source Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << efield3 << endl;
	      cout << "Source term is:                   " << efield2 << endl;
	      cout << "MASA term is:                     " << efield << endl;
	      cout << x << " " << y << endl;
	      exit(1);
	    }

	  if(exact_p3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "P Field Analytical Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << exact_p << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(rho3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "RHO Source Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << rho << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(exact_rho3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "RHO Analytical Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << exact_rho << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  // adding a new error check: ensure physical results are coming out!
	  if(0 > exact_rho)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "Initial Variables are returning non-physical results!\n";
	      cout << "RHO analytical\n";
	      exit(1);
	    }
	  
	  if(0 > exact_p)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "Initial Variables are returning non-physical results!\n";
	      cout << "Pressure is negative!\n";
	    exit(1);
	    }
	  
	}// done iterating
  
  // reroute stdout for regressions: TODO remove when logger mechanism
  // is used inside masa; these tests currently just verify functions
  // run successfully.
  freopen("/dev/null","w",stdout);

  // test gradient error terms
  Scalar derr = masa_eval_grad_u<Scalar>(0,0,0,0);
  if(derr != -1)
    {
      cout << "MASA :: gradient (0) error condition failed!\n";
      exit(1);
    }
  
  derr = masa_eval_grad_u<Scalar>(0,0,0,4);
  if(derr != -1)
    {
      cout << "MASA :: gradient (4) error condition failed!\n";
      exit(1);
    }

  // v
  derr = masa_eval_grad_v<Scalar>(0,0,0,0);
  if(derr != -1)
    {
      cout << "MASA :: v gradient (0) error condition failed!\n";
      exit(1);
    }
  
  derr = masa_eval_grad_v<Scalar>(0,0,0,4);
  if(derr != -1)
    {
      cout << "MASA :: v gradient (4) error condition failed!\n";
      exit(1);
    }

  // w
  derr = masa_eval_grad_w<Scalar>(0,0,0,0);
  if(derr != -1)
    {
      cout << "MASA :: w gradient (0) error condition failed!\n";
      exit(1);
    }
  
  derr = masa_eval_grad_w<Scalar>(0,0,0,4);
  if(derr != -1)
    {
      cout << "MASA :: w gradient (4) error condition failed!\n";
      exit(1);
    }

  // p
  derr = masa_eval_grad_p<Scalar>(0,0,0,0);
  if(derr != -1)
    {
      cout << "MASA :: p gradient (0) error condition failed!\n";
      exit(1);
    }
  
  derr = masa_eval_grad_p<Scalar>(0,0,0,4);
  if(derr != -1)
    {
      cout << "MASA :: p gradient (4) error condition failed!\n";
      exit(1);
    }

  // rho
  derr = masa_eval_grad_rho<Scalar>(0,0,0,0);
  if(derr != -1)
    {
      cout << "MASA :: rho gradient (0) error condition failed!\n";
      exit(1);
    }
  
  derr = masa_eval_grad_rho<Scalar>(0,0,0,4);
  if(derr != -1)
    {
      cout << "MASA :: rho gradient (4) error condition failed!\n";
      exit(1);
    }

  // all tests passed
  return 0;
}

int main()
{
  int err=0;

  err += run_regression<double>();
  //err += run_regression<long double>();

  return err;
}
