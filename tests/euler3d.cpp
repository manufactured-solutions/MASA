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

#include <tests.h>
#include <cmath>

using namespace MASA;
using namespace std;

typedef double Scalar;

template<typename Scalar>
Scalar anQ_p (Scalar x,Scalar y,Scalar z,Scalar p_0,Scalar p_x,Scalar p_y,Scalar p_z,Scalar a_px,Scalar a_py,Scalar a_pz,Scalar L)
{
  Scalar pi = std::acos(Scalar(-1));
  Scalar exact_p = p_0 + p_x * std::cos(a_px * pi * x / L) + p_y * std::sin(a_py * pi * y / L) + p_z * std::cos(a_pz * pi * z / L);
  return exact_p;
}
 
template<typename Scalar>
Scalar anQ_u (Scalar x,Scalar y,Scalar z,Scalar u_0,Scalar u_x,Scalar u_y,Scalar u_z,Scalar a_ux,Scalar a_uy,Scalar a_uz,Scalar L)
{
  Scalar pi = std::acos(Scalar(-1));
  Scalar exact_u = u_0 + u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L) + u_z * std::cos(a_uz * pi * z / L);  
  return exact_u;
} 

template<typename Scalar>
Scalar anQ_v (Scalar x,Scalar y,Scalar z,Scalar v_0,Scalar v_x,Scalar v_y,Scalar v_z,Scalar a_vx,Scalar a_vy,Scalar a_vz,Scalar L)
{
  Scalar pi = std::acos(Scalar(-1));
  Scalar exact_v = v_0 + v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L) + v_z * std::sin(a_vz * pi * z / L);
  return exact_v;
}

template<typename Scalar>
Scalar anQ_w (Scalar x,Scalar y,Scalar z,Scalar w_0,Scalar w_x,Scalar w_y,Scalar w_z,Scalar a_wx,Scalar a_wy,Scalar a_wz,Scalar L)
{
  Scalar pi = std::acos(Scalar(-1));
  Scalar exact_w = w_0 + w_x * std::sin(a_wx * pi * x / L) + w_y * std::sin(a_wy * pi * y / L) + w_z * std::cos(a_wz * pi * z / L);
  return exact_w;
}

template<typename Scalar>
Scalar anQ_rho (Scalar x,Scalar y,Scalar z,Scalar rho_0,Scalar rho_x,Scalar rho_y,Scalar rho_z,Scalar a_rhox,Scalar a_rhoy,Scalar a_rhoz,Scalar L)
{ 
  Scalar pi = std::acos(Scalar(-1));
  Scalar exact_rho = rho_0 + rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L) + rho_z * std::sin(a_rhoz * pi * z / L);
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
  Scalar pi = std::acos(Scalar(-1));
  Scalar Q_e;
  Scalar RHO;
  Scalar P;
  Scalar U;
  Scalar V;
  Scalar W;

  RHO = rho_0 + rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L) + rho_z * std::sin(a_rhoz * pi * z / L);
  U = u_0 + u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L) + u_z * std::cos(a_uz * pi * z / L);
  V = v_0 + v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L) + v_z * std::sin(a_vz * pi * z / L);
  W = w_0 + w_x * std::sin(a_wx * pi * x / L) + w_y * std::sin(a_wy * pi * y / L) + w_z * std::cos(a_wz * pi * z / L);
  P = p_0 + p_x * std::cos(a_px * pi * x / L) + p_y * std::sin(a_py * pi * y / L) + p_z * std::cos(a_pz * pi * z / L);

  Q_e = -a_px * pi * p_x * Gamma * U * std::sin(a_px * pi * x / L) / (Gamma - Scalar(0.1e1)) / L + a_py * pi * p_y * Gamma * V * std::cos(a_py * pi * y / L) / (Gamma - Scalar(0.1e1)) / L - a_pz * pi * p_z * Gamma * W * std::sin(a_pz * pi * z / L) / (Gamma - Scalar(0.1e1)) / L + (U * U + V * V + W * W) * a_rhox * pi * rho_x * U * std::cos(a_rhox * pi * x / L) / L / Scalar(0.2e1) - (U * U + V * V + W * W) * a_rhoy * pi * rho_y * V * std::sin(a_rhoy * pi * y / L) / L / Scalar(0.2e1) + (U * U + V * V + W * W) * a_rhoz * pi * rho_z * W * std::cos(a_rhoz * pi * z / L) / L / Scalar(0.2e1) - (Scalar(-0.3e1) * a_ux * u_x * std::cos(a_ux * pi * x / L) - a_vy * v_y * std::cos(a_vy * pi * y / L) + a_wz * w_z * std::sin(a_wz * pi * z / L)) * pi * RHO * U * U / L / Scalar(0.2e1) - (a_uy * u_y * std::sin(a_uy * pi * y / L) + a_vx * v_x * std::sin(a_vx * pi * x / L)) * pi * RHO * U * V / L - (a_uz * u_z * std::sin(a_uz * pi * z / L) - a_wx * w_x * std::cos(a_wx * pi * x / L)) * pi * RHO * U * W / L - (-a_ux * u_x * std::cos(a_ux * pi * x / L) - Scalar(0.3e1) * a_vy * v_y * std::cos(a_vy * pi * y / L) + a_wz * w_z * std::sin(a_wz * pi * z / L)) * pi * RHO * V * V / L / Scalar(0.2e1) + (a_vz * v_z * std::cos(a_vz * pi * z / L) + a_wy * w_y * std::cos(a_wy * pi * y / L)) * pi * RHO * V * W / L - (-a_ux * u_x * std::cos(a_ux * pi * x / L) - a_vy * v_y * std::cos(a_vy * pi * y / L) + Scalar(0.3e1) * a_wz * w_z * std::sin(a_wz * pi * z / L)) * pi * RHO * W * W / L / Scalar(0.2e1) - (-a_ux * u_x * std::cos(a_ux * pi * x / L) - a_vy * v_y * std::cos(a_vy * pi * y / L) + a_wz * w_z * std::sin(a_wz * pi * z / L)) * pi * Gamma * P / (Gamma - Scalar(0.1e1)) / L;

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
  Scalar p_x,
  Scalar a_px,
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
  Scalar pi = std::acos(Scalar(-1));
  Scalar Q_u;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar W;

  RHO = rho_0 + rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L) + rho_z * std::sin(a_rhoz * pi * z / L);
  U = u_0 + u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L) + u_z * std::cos(a_uz * pi * z / L);
  V = v_0 + v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L) + v_z * std::sin(a_vz * pi * z / L);
  W = w_0 + w_x * std::sin(a_wx * pi * x / L) + w_y * std::sin(a_wy * pi * y / L) + w_z * std::cos(a_wz * pi * z / L);

  Q_u = a_rhox * pi * rho_x * U * U * std::cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * U * V * std::sin(a_rhoy * pi * y / L) / L + a_rhoz * pi * rho_z * U * W * std::cos(a_rhoz * pi * z / L) / L - a_uy * pi * u_y * RHO * V * std::sin(a_uy * pi * y / L) / L - a_uz * pi * u_z * RHO * W * std::sin(a_uz * pi * z / L) / L - a_px * pi * p_x * std::sin(a_px * pi * x / L) / L + (Scalar(0.2e1) * a_ux * u_x * std::cos(a_ux * pi * x / L) + a_vy * v_y * std::cos(a_vy * pi * y / L) - a_wz * w_z * std::sin(a_wz * pi * z / L)) * pi * RHO * U / L;

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
  Scalar p_y,
  Scalar p_z,
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
  (void) p_z;
  (void) a_pz;

  Scalar pi = std::acos(Scalar(-1));
  Scalar Q_v;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar W;

  RHO = rho_0 + rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L) + rho_z * std::sin(a_rhoz * pi * z / L);
  U = u_0 + u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L) + u_z * std::cos(a_uz * pi * z / L);
  V = v_0 + v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L) + v_z * std::sin(a_vz * pi * z / L);
  W = w_0 + w_x * std::sin(a_wx * pi * x / L) + w_y * std::sin(a_wy * pi * y / L) + w_z * std::cos(a_wz * pi * z / L);

  Q_v = a_rhox * pi * rho_x * U * V * std::cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * V * V * std::sin(a_rhoy * pi * y / L) / L + a_rhoz * pi * rho_z * V * W * std::cos(a_rhoz * pi * z / L) / L - a_vx * pi * v_x * RHO * U * std::sin(a_vx * pi * x / L) / L + a_vz * pi * v_z * RHO * W * std::cos(a_vz * pi * z / L) / L + a_py * pi * p_y * std::cos(a_py * pi * y / L) / L + (a_ux * u_x * std::cos(a_ux * pi * x / L) + Scalar(0.2e1) * a_vy * v_y * std::cos(a_vy * pi * y / L) - a_wz * w_z * std::sin(a_wz * pi * z / L)) * pi * RHO * V / L;

  return(Q_v);
}

template<typename Scalar>
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
  Scalar p_z,
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
  Scalar pi = std::acos(Scalar(-1));
  Scalar Q_w;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar W;

  RHO = rho_0 + rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L) + rho_z * std::sin(a_rhoz * pi * z / L);
  U = u_0 + u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L) + u_z * std::cos(a_uz * pi * z / L);
  V = v_0 + v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L) + v_z * std::sin(a_vz * pi * z / L);
  W = w_0 + w_x * std::sin(a_wx * pi * x / L) + w_y * std::sin(a_wy * pi * y / L) + w_z * std::cos(a_wz * pi * z / L);

  Q_w = a_rhox * pi * rho_x * U * W * std::cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * V * W * std::sin(a_rhoy * pi * y / L) / L + a_rhoz * pi * rho_z * W * W * std::cos(a_rhoz * pi * z / L) / L + a_wx * pi * w_x * RHO * U * std::cos(a_wx * pi * x / L) / L + a_wy * pi * w_y * RHO * V * std::cos(a_wy * pi * y / L) / L - a_pz * pi * p_z * std::sin(a_pz * pi * z / L) / L + (a_ux * u_x * std::cos(a_ux * pi * x / L) + a_vy * v_y * std::cos(a_vy * pi * y / L) - Scalar(0.2e1) * a_wz * w_z * std::sin(a_wz * pi * z / L)) * pi * RHO * W / L;

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
  Scalar p_y,
  Scalar p_z,
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
  Scalar pi = std::acos(Scalar(-1));
  Scalar Q_rho;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar W;
  (void) p_y;
  (void) a_py;
  (void) p_z;
  (void) a_pz;
  (void) mu;


  RHO = rho_0 + rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L) + rho_z * std::sin(a_rhoz * pi * z / L);
  U = u_0 + u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L) + u_z * std::cos(a_uz * pi * z / L);
  V = v_0 + v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L) + v_z * std::sin(a_vz * pi * z / L);
  W = w_0 + w_x * std::sin(a_wx * pi * x / L) + w_y * std::sin(a_wy * pi * y / L) + w_z * std::cos(a_wz * pi * z / L);

  Q_rho = a_rhox * pi * rho_x * U * std::cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * V * std::sin(a_rhoy * pi * y / L) / L + a_rhoz * pi * rho_z * W * std::cos(a_rhoz * pi * z / L) / L + (a_ux * u_x * std::cos(a_ux * pi * x / L) + a_vy * v_y * std::cos(a_vy * pi * y / L) - a_wz * w_z * std::sin(a_wz * pi * z / L)) * pi * RHO / L;

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
	  efield = masa_eval_source_rho_e  <Scalar>(x,y,z);
	  rho    = masa_eval_source_rho    <Scalar>(x,y,z);
	  
	  //evaluate analytical terms
	  exact_u = masa_eval_exact_u      <Scalar>(x,y,z);
	  exact_v = masa_eval_exact_v      <Scalar>(x,y,z);
	  exact_w = masa_eval_exact_w      <Scalar>(x,y,z);
	  exact_p = masa_eval_exact_p      <Scalar>(x,y,z);
	  exact_rho = masa_eval_exact_rho  <Scalar>(x,y,z);

	  // eval gradient terms
	  gradx = masa_eval_grad_u<Scalar>(x,y,z,1);
	  nancheck(gradx);
	  grady = masa_eval_grad_u<Scalar>(x,y,z,2);
	  nancheck(grady);	
	  gradz = masa_eval_grad_u<Scalar>(x,y,z,3);
	  nancheck(gradz);

	  gradx = masa_eval_grad_v<Scalar>(x,y,z,1);
	  nancheck(gradx);
	  grady = masa_eval_grad_v<Scalar>(x,y,z,2);
	  nancheck(grady);	
	  gradz = masa_eval_grad_v<Scalar>(x,y,z,3);
	  nancheck(gradz);

	  gradx = masa_eval_grad_w<Scalar>(x,y,z,1);
	  nancheck(gradx);
	  grady = masa_eval_grad_w<Scalar>(x,y,z,2);
	  nancheck(grady);	
	  gradz = masa_eval_grad_w<Scalar>(x,y,z,3);
	  nancheck(gradz);

	  gradx = masa_eval_grad_p<Scalar>(x,y,z,1);
	  nancheck(gradx);
	  grady = masa_eval_grad_p<Scalar>(x,y,z,2);
	  nancheck(grady);	
	  gradz = masa_eval_grad_p<Scalar>(x,y,z,3);
	  nancheck(gradz);

	  gradx = masa_eval_grad_rho<Scalar>(x,y,z,1);
	  nancheck(gradx);
	  grady = masa_eval_grad_rho<Scalar>(x,y,z,2);
	  nancheck(grady);	
	  gradz = masa_eval_grad_rho<Scalar>(x,y,z,3);
	  nancheck(gradz);

	  // check against maple output
	  ufield2   = SourceQ_u  (x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_x,a_px,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,L);
	  vfield2   = SourceQ_v  (x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_y,p_z,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,L);
	  wfield2   = SourceQ_w  (x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_z,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,L);
	  rho2      = SourceQ_rho(x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_y,p_z,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,mu,L);
	  efield2   = SourceQ_e  (x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,mu,Gamma,L);
  
	  exact_u2     = anQ_u   (x,y,z,u_0,u_x,u_y,u_z,a_ux,a_uy,a_uz,L);
	  exact_v2     = anQ_v   (x,y,z,v_0,v_x,v_y,v_z,a_vx,a_vy,a_vz,L);
	  exact_w2     = anQ_w   (x,y,z,w_0,w_x,w_y,w_z,a_wx,a_wy,a_wz,L);
	  exact_rho2   = anQ_rho (x,y,z,rho_0,rho_x,rho_y,rho_z,a_rhox,a_rhoy,a_rhoz,L);
	  exact_p2     = anQ_p   (x,y,z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,L);

	  // test the result is roughly zero
	  // choose between abs and rel error
#ifdef MASA_STRICT_REGRESSION

	  ufield3 = std::abs(ufield-ufield2);
	  vfield3 = std::abs(vfield-vfield2);
	  wfield3 = std::abs(wfield-wfield2);
	  efield3 = std::abs(efield-efield2);
	  rho3    = std::abs(rho-rho2);

	  exact_u3   = std::abs(exact_u-exact_u2);
	  exact_v3   = std::abs(exact_v-exact_v2);
	  exact_w3   = std::abs(exact_w-exact_w2);
	  exact_rho3 = std::abs(exact_rho-exact_rho2);
	  exact_p3   = std::abs(exact_p-exact_p2);

#else

	  ufield3 = std::abs(ufield-ufield2)/std::abs(ufield2);
	  vfield3 = std::abs(vfield-vfield2)/std::abs(vfield2);
	  wfield3 = std::abs(wfield-wfield2)/std::abs(wfield2);
	  efield3 = std::abs(efield-efield2)/std::abs(efield2);
	  rho3    = std::abs(rho-rho2)/std::abs(rho2);

	  exact_u3   = std::abs(exact_u-exact_u2)/std::abs(exact_u2);
	  exact_v3   = std::abs(exact_v-exact_v2)/std::abs(exact_v2);
	  exact_w3   = std::abs(exact_w-exact_w2)/std::abs(exact_w2);
	  exact_rho3 = std::abs(exact_rho-exact_rho2)/std::abs(exact_rho2);
	  exact_p3   = std::abs(exact_p-exact_p2)/std::abs(exact_p2);
#endif	  

	  threshcheck(ufield3,threshold);
	  threshcheck(vfield3,threshold);
	  threshcheck(wfield3,threshold);
	  threshcheck(efield3,threshold);
	  threshcheck(rho3,threshold);
	  
	  threshcheck(exact_u3,threshold);
	  threshcheck(exact_v3,threshold);
	  threshcheck(exact_w3,threshold);
	  threshcheck(exact_rho3,threshold);
	  threshcheck(exact_p3,threshold);

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
  test_grad(derr);

  derr = masa_eval_grad_u<Scalar>(0,0,0,4);
  test_grad(derr);

  // v
  derr = masa_eval_grad_v<Scalar>(0,0,0,0);
  test_grad(derr);

  derr = masa_eval_grad_v<Scalar>(0,0,0,4);
  test_grad(derr);

  // w
  derr = masa_eval_grad_w<Scalar>(0,0,0,0);
  test_grad(derr);
  
  derr = masa_eval_grad_w<Scalar>(0,0,0,4);
  test_grad(derr);

  // p
  derr = masa_eval_grad_p<Scalar>(0,0,0,0);
  test_grad(derr);
  
  derr = masa_eval_grad_p<Scalar>(0,0,0,4);
  test_grad(derr);

  // rho
  derr = masa_eval_grad_rho<Scalar>(0,0,0,0);
  test_grad(derr);
  
  derr = masa_eval_grad_rho<Scalar>(0,0,0,4);
  test_grad(derr);

  // all tests passed
  return 0;
}

int main()
{
  int err=0;

  err += run_regression<double>();
  err += run_regression<long double>();

  return err;
}
