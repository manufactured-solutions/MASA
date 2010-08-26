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

const double PI = acos(-1);

double SourceQ_u (
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
  double Q_u;
  Q_u = 0.4e1 / 0.3e1 * mu * u_x * sin(a_ux * PI * x / L) * a_ux * a_ux * PI * PI * pow(L, -0.2e1) + mu * u_y * cos(a_uy * PI * y / L) * a_uy * a_uy * PI * PI * pow(L, -0.2e1) + mu * u_z * cos(a_uz * PI * z / L) * a_uz * a_uz * PI * PI * pow(L, -0.2e1) - p_x * sin(a_px * PI * x / L) * a_px * PI / L + rho_x * cos(a_rhox * PI * x / L) * pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_rhoy * PI / L + rho_z * cos(a_rhoz * PI * z / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_rhoz * PI / L + 0.2e1 * u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_ux * PI / L - u_y * sin(a_uy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_uy * PI / L - u_z * sin(a_uz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_uz * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_vy * PI / L - w_z * sin(a_wz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_wz * PI / L;
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
  double mu,
  double L,
  double R,
  double k)
{
  double Q_v;
  Q_v = mu * v_x * cos(a_vx * PI * x / L) * a_vx * a_vx * PI * PI * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * mu * v_y * sin(a_vy * PI * y / L) * a_vy * a_vy * PI * PI * pow(L, -0.2e1) + mu * v_z * sin(a_vz * PI * z / L) * a_vz * a_vz * PI * PI * pow(L, -0.2e1) + p_y * cos(a_py * PI * y / L) * a_py * PI / L + rho_x * cos(a_rhox * PI * x / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L), 0.2e1) * a_rhoy * PI / L + rho_z * cos(a_rhoz * PI * z / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_rhoz * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_ux * PI / L - v_x * sin(a_vx * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_vx * PI / L + 0.2e1 * v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_vy * PI / L + v_z * cos(a_vz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_vz * PI / L - w_z * sin(a_wz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_wz * PI / L;
  return(Q_v);
}

#include <math.h>

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
  double mu,
  double L,
  double R,
  double k)
{
  double Q_w;
  Q_w = mu * w_x * sin(a_wx * PI * x / L) * a_wx * a_wx * PI * PI * pow(L, -0.2e1) + mu * w_y * sin(a_wy * PI * y / L) * a_wy * a_wy * PI * PI * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * mu * w_z * cos(a_wz * PI * z / L) * a_wz * a_wz * PI * PI * pow(L, -0.2e1) - p_z * sin(a_pz * PI * z / L) * a_pz * PI / L + rho_x * cos(a_rhox * PI * x / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_rhoy * PI / L + rho_z * cos(a_rhoz * PI * z / L) * pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1) * a_rhoz * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_ux * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_vy * PI / L + w_x * cos(a_wx * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_wx * PI / L + w_y * cos(a_wy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_wy * PI / L - 0.2e1 * w_z * sin(a_wz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_wz * PI / L;
  return(Q_w);
}


double SourceQ_e (
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
    
    // correct version
    const double Q_e = -sin(a_rhoy * PI * y / L) * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0) * (pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1) + pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0, 0.2e1)) * rho_y * a_rhoy * PI / L / 0.2e1 + cos(a_rhox * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * (pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1) + pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0, 0.2e1)) * rho_x * a_rhox * PI / L / 0.2e1 + 0.4e1 / 0.3e1 * (-pow(cos(a_ux * PI * x / L), 0.2e1) * u_x + sin(a_ux * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L))) * mu * u_x * a_ux * a_ux * PI * PI * pow(L, -0.2e1) + (-pow(sin(a_uy * PI * y / L), 0.2e1) * u_y + cos(a_uy * PI * y / L) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L))) * mu * u_y * a_uy * a_uy * PI * PI * pow(L, -0.2e1) + (-pow(sin(a_uz * PI * z / L), 0.2e1) * u_z + cos(a_uz * PI * z / L) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L))) * mu * u_z * a_uz * a_uz * PI * PI * pow(L, -0.2e1) - (pow(sin(a_vx * PI * x / L), 0.2e1) * v_x - cos(a_vx * PI * x / L) * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0)) * mu * v_x * a_vx * a_vx * PI * PI * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * (pow(cos(a_vy * PI * y / L), 0.2e1) * v_y - sin(a_vy * PI * y / L) * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0)) * mu * v_y * a_vy * a_vy * PI * PI * pow(L, -0.2e1) - (pow(cos(a_vz * PI * z / L), 0.2e1) * v_z - sin(a_vz * PI * z / L) * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0)) * mu * v_z * a_vz * a_vz * PI * PI * pow(L, -0.2e1) + (-pow(cos(a_wx * PI * x / L), 0.2e1) * w_x + sin(a_wx * PI * x / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L))) * mu * w_x * a_wx * a_wx * PI * PI * pow(L, -0.2e1) + (-pow(cos(a_wy * PI * y / L), 0.2e1) * w_y + sin(a_wy * PI * y / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L))) * mu * w_y * a_wy * a_wy * PI * PI * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * (-pow(sin(a_wz * PI * z / L), 0.2e1) * w_z + cos(a_wz * PI * z / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L))) * mu * w_z * a_wz * a_wz * PI * PI * pow(L, -0.2e1) + (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0) * w_y * cos(a_wy * PI * y / L) * a_wy * PI * (-(0.3e1 * pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1) + pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0, 0.2e1) + pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1)) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) / L / 0.2e1 - Gamma * (p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L)) / L / (Gamma - 0.1e1)) * w_z * sin(a_wz * PI * z / L) * a_wz * PI / L + (-0.1e1) * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * v_x * sin(a_vx * PI * x / L) * a_vx * PI * ((pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1) + 0.3e1 * pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0, 0.2e1)) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) / L / 0.2e1 + Gamma * (p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L)) / L / (Gamma - 0.1e1)) * v_y * cos(a_vy * PI * y / L) * a_vy * PI / L - 0.2e1 * cos(a_rhox * PI * x / L) * rho_x * sin(a_px * PI * x / L) * k * p_x * a_px * a_rhox * PI * PI * pow(L, -0.2e1) / R * pow(rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L), -0.2e1) + sin(a_py * PI * y / L) * k * p_y * a_py * a_py * PI * PI * pow(L, -0.2e1) / R / (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) - 0.2e1 * sin(a_rhoy * PI * y / L) * rho_y * cos(a_py * PI * y / L) * k * p_y * a_py * a_rhoy * PI * PI * pow(L, -0.2e1) / R * pow(rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L), -0.2e1) - 0.2e1 * mu * v_z * w_y * cos(a_vz * PI * z / L) * cos(a_wy * PI * y / L) * a_vz * a_wy * PI * PI * pow(L, -0.2e1) + cos(a_pz * PI * z / L) * k * p_z * a_pz * a_pz * PI * PI * pow(L, -0.2e1) / R / (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) + cos(a_px * PI * x / L) * k * p_x * a_px * a_px * PI * PI * pow(L, -0.2e1) / R / (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) + (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0) * v_z * cos(a_vz * PI * z / L) * a_vz * PI / L - (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * u_y * sin(a_uy * PI * y / L) * a_uy * PI / L + (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * w_x * cos(a_wx * PI * x / L) * a_wx * PI / L - (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * u_z * sin(a_uz * PI * z / L) * a_uz * PI / L - 0.2e1 * cos(a_rhoz * PI * z / L) * rho_z * sin(a_pz * PI * z / L) * k * p_z * a_pz * a_rhoz * PI * PI * pow(L, -0.2e1) / R * pow(rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L), -0.2e1) - (0.2e1 * pow(cos(a_rhox * PI * x / L), 0.2e1) * rho_x + sin(a_rhox * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L))) * (p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L)) * k * rho_x * a_rhox * a_rhox * PI * PI * pow(L, -0.2e1) / R * pow(rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L), -0.3e1) - (0.2e1 * pow(sin(a_rhoy * PI * y / L), 0.2e1) * rho_y + cos(a_rhoy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L))) * (p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L)) * k * rho_y * a_rhoy * a_rhoy * PI * PI * pow(L, -0.2e1) / R * pow(rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L), -0.3e1) - (0.2e1 * pow(cos(a_rhoz * PI * z / L), 0.2e1) * rho_z + sin(a_rhoz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L))) * (p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L)) * k * rho_z * a_rhoz * a_rhoz * PI * PI * pow(L, -0.2e1) / R * pow(rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L), -0.3e1) + 0.4e1 / 0.3e1 * mu * u_x * v_y * cos(a_ux * PI * x / L) * cos(a_vy * PI * y / L) * a_ux * a_vy * PI * PI * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * mu * u_x * w_z * cos(a_ux * PI * x / L) * sin(a_wz * PI * z / L) * a_ux * a_wz * PI * PI * pow(L, -0.2e1) - 0.2e1 * mu * u_y * v_x * sin(a_uy * PI * y / L) * sin(a_vx * PI * x / L) * a_uy * a_vx * PI * PI * pow(L, -0.2e1) + 0.2e1 * mu * u_z * w_x * cos(a_wx * PI * x / L) * sin(a_uz * PI * z / L) * a_uz * a_wx * PI * PI * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * mu * v_y * w_z * cos(a_vy * PI * y / L) * sin(a_wz * PI * z / L) * a_vy * a_wz * PI * PI * pow(L, -0.2e1) + 0.1e1 / 0.2e1 * cos(a_rhoz * PI * z / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1) + pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0, 0.2e1)) * rho_z * a_rhoz * PI * ((pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1) + pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0, 0.2e1) + 0.3e1 * pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1)) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) / L / 0.2e1 + Gamma * (p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L)) / L / (Gamma - 0.1e1)) * u_x * cos(a_ux * PI * x / L) * a_ux * PI / L - Gamma * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * sin(a_px * PI * x / L) * p_x * a_px * PI / L / (Gamma - 0.1e1) + Gamma * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0) * cos(a_py * PI * y / L) * p_y * a_py * PI / L / (Gamma - 0.1e1) - Gamma * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * sin(a_pz * PI * z / L) * p_z * a_pz * PI / L / (Gamma - 0.1e1);
    
    // juans version 
    /*
    const double unique2 = pow(L,(-1));
    const double unique63 = (a_px * PI * x * unique2);
    const double d_p_d_x = (-a_px * p_x * PI * unique2 * sin(unique63));
    const double unique1 = (a_rhox * PI * x * unique2);
    const double d_rho_d_x = (a_rhox * PI * rho_x * unique2 * cos(unique1));
    const double unique62 = cos(unique63);
    const double unique65 = (a_py * PI * y * unique2);
    const double unique64 = sin(unique65);
    const double unique67 = (a_pz * PI * z * unique2);
    const double unique66 = cos(unique67);
    const double p = (p_0 + (p_x * unique62) + (p_y * unique64) + (p_z * unique66));
    const double unique0 = sin(unique1);
    const double unique41 = (a_rhoy * PI * y * unique2);
    const double unique40 = cos(unique41);
    const double unique43 = (a_rhoz * PI * z * unique2);
    const double unique42 = sin(unique43);
    const double rho = (rho_0 + (rho_x * unique0) + (rho_y * unique40) + (rho_z * unique42));
    const double unique32 = (R * rho);
    const double unique31 = pow(unique32,(-1));
    const double unique68 = pow(unique32,(-2));
    const double d_T_d_x = (-(R * d_rho_d_x * p * unique68) + (d_p_d_x * unique31));
    const double unique34 = (-1 + Gamma);
    const double unique33 = pow(unique34,(-1));
    const double d_e_d_x = (R * d_T_d_x * unique33);
    const double unique45 = (a_ux * PI * x * unique2);
    const double d_u_d_x = (a_ux * PI * u_x * unique2 * cos(unique45));
    const double unique51 = (a_vx * PI * x * unique2);
    const double d_v_d_x = (-a_vx * PI * v_x * unique2 * sin(unique51));
    const double unique57 = (a_wx * PI * x * unique2);
    const double d_w_d_x = (a_wx * PI * w_x * unique2 * cos(unique57));
    const double unique44 = sin(unique45);
    const double unique47 = (a_uy * PI * y * unique2);
    const double unique46 = cos(unique47);
    const double unique49 = (a_uz * PI * z * unique2);
    const double unique48 = cos(unique49);
    const double u = (u_0 + (u_x * unique44) + (u_y * unique46) + (u_z * unique48));
    const double unique50 = cos(unique51);
    const double unique53 = (a_vy * PI * y * unique2);
    const double unique52 = sin(unique53);
    const double unique55 = (a_vz * PI * z * unique2);
    const double unique54 = sin(unique55);
    const double v = (v_0 + (v_x * unique50) + (v_y * unique52) + (v_z * unique54));
    const double unique56 = sin(unique57);
    const double unique59 = (a_wy * PI * y * unique2);
    const double unique58 = sin(unique59);
    const double unique61 = (a_wz * PI * z * unique2);
    const double unique60 = cos(unique61);
    const double w = (w_0 + (w_x * unique56) + (w_y * unique58) + (w_z * unique60));
    const double d_e_t_d_x = ((5.000000000000000e-01 * ((2 * d_u_d_x * u) + (2 * d_v_d_x * v) + (2 * d_w_d_x * w))) + d_e_d_x);
    const double d_p_d_y = (a_py * p_y * PI * unique2 * cos(unique65));
    const double d_rho_d_y = (-a_rhoy * PI * rho_y * unique2 * sin(unique41));
    const double d_T_d_y = (-(R * d_rho_d_y * p * unique68) + (d_p_d_y * unique31));
    const double d_e_d_y = (R * d_T_d_y * unique33);
    const double d_u_d_y = (-a_uy * PI * u_y * unique2 * sin(unique47));
    const double d_v_d_y = (a_vy * PI * v_y * unique2 * cos(unique53));
    const double d_w_d_y = (a_wy * PI * w_y * unique2 * cos(unique59));
    const double d_e_t_d_y = ((5.000000000000000e-01 * ((2 * d_u_d_y * u) + (2 * d_v_d_y * v) + (2 * d_w_d_y * w))) + d_e_d_y);
    const double d_p_d_z = (-a_pz * p_z * PI * unique2 * sin(unique67));
    const double d_rho_d_z = (a_rhoz * PI * rho_z * unique2 * cos(unique43));
    const double d_T_d_z = (-(R * d_rho_d_z * p * unique68) + (d_p_d_z * unique31));
    const double d_e_d_z = (R * d_T_d_z * unique33);
    const double d_u_d_z = (-a_uz * PI * u_z * unique2 * sin(unique49));
    const double d_v_d_z = (a_vz * PI * v_z * unique2 * cos(unique55));
    const double d_w_d_z = (-a_wz * PI * w_z * unique2 * sin(unique61));
    const double d_e_t_d_z = ((5.000000000000000e-01 * ((2 * d_u_d_z * u) + (2 * d_v_d_z * v) + (2 * d_w_d_z * w))) + d_e_d_z);
    const double unique36 = pow(L,(-2));
    const double unique37 = pow(PI,2);
    const double d_d_p_d_x_d_x = (-p_x * pow(a_px,2) * unique36 * unique37 * unique62);
    const double d_d_rho_d_x_d_x = (-rho_x * pow(a_rhox,2) * unique0 * unique36 * unique37);
    const double unique39 = pow(R,2);
    const double unique69 = pow(unique32,(-3));
    const double d_d_T_d_x_d_x = (-(2 * R * d_p_d_x * d_rho_d_x * unique68) - (R * d_d_rho_d_x_d_x * p * unique68) + (2 * pow(d_rho_d_x,2) * p * unique39 * unique69) + (d_d_p_d_x_d_x * unique31));
    const double d_q_x_d_x = (-k * d_d_T_d_x_d_x);
    const double d_d_p_d_y_d_y = (-p_y * pow(a_py,2) * unique36 * unique37 * unique64);
    const double d_d_rho_d_y_d_y = (-rho_y * pow(a_rhoy,2) * unique36 * unique37 * unique40);
    const double d_d_T_d_y_d_y = (-(2 * R * d_p_d_y * d_rho_d_y * unique68) - (R * d_d_rho_d_y_d_y * p * unique68) + (2 * pow(d_rho_d_y,2) * p * unique39 * unique69) + (d_d_p_d_y_d_y * unique31));
    const double d_q_y_d_y = (-k * d_d_T_d_y_d_y);
    const double d_d_p_d_z_d_z = (-p_z * pow(a_pz,2) * unique36 * unique37 * unique66);
    const double d_d_rho_d_z_d_z = (-rho_z * pow(a_rhoz,2) * unique36 * unique37 * unique42);
    const double d_d_T_d_z_d_z = (-(2 * R * d_p_d_z * d_rho_d_z * unique68) - (R * d_d_rho_d_z_d_z * p * unique68) + (2 * pow(d_rho_d_z,2) * p * unique39 * unique69) + (d_d_p_d_z_d_z * unique31));
    const double d_q_z_d_z = (-k * d_d_T_d_z_d_z);
    const double d_d_u_d_x_d_x = (-u_x * pow(a_ux,2) * unique36 * unique37 * unique44);
    const double d_tau_xx_d_x = (1.333333333333333e+00 * mu * d_d_u_d_x_d_x);
    const double d_d_v_d_x_d_x = (-v_x * pow(a_vx,2) * unique36 * unique37 * unique50);
    const double d_tau_xy_d_x = (mu * d_d_v_d_x_d_x);
    const double d_d_w_d_x_d_x = (-w_x * pow(a_wx,2) * unique36 * unique37 * unique56);
    const double d_tau_xz_d_x = (mu * d_d_w_d_x_d_x);
    const double d_d_u_d_y_d_y = (-u_y * pow(a_uy,2) * unique36 * unique37 * unique46);
    const double d_tau_xy_d_y = (mu * d_d_u_d_y_d_y);
    const double d_tau_yx_d_y = d_tau_xy_d_y;
    const double d_d_v_d_y_d_y = (-v_y * pow(a_vy,2) * unique36 * unique37 * unique52);
    const double d_tau_yy_d_y = (1.333333333333333e+00 * mu * d_d_v_d_y_d_y);
    const double d_d_w_d_y_d_y = (-w_y * pow(a_wy,2) * unique36 * unique37 * unique58);
    const double d_tau_yz_d_y = (mu * d_d_w_d_y_d_y);
    const double d_d_u_d_z_d_z = (-u_z * pow(a_uz,2) * unique36 * unique37 * unique48);
    const double d_tau_xz_d_z = (mu * d_d_u_d_z_d_z);
    const double d_tau_zx_d_z = d_tau_xz_d_z;
    const double d_d_v_d_z_d_z = (-v_z * pow(a_vz,2) * unique36 * unique37 * unique54);
    const double d_tau_yz_d_z = (mu * d_d_v_d_z_d_z);
    const double d_tau_zy_d_z = d_tau_yz_d_z;
    const double d_d_w_d_z_d_z = (-w_z * pow(a_wz,2) * unique36 * unique37 * unique60);
    const double d_tau_zz_d_z = (1.333333333333333e+00 * mu * d_d_w_d_z_d_z);
    const double T = (p * unique31);
    const double e = (R * T * unique33);
    const double e_t = ((5.000000000000000e-01 * (pow(u,2) + pow(v,2) + pow(w,2))) + e);
    const double tau_xx = (6.666666666666666e-01 * mu * (-d_v_d_y - d_w_d_z + (2 * d_u_d_x)));
    const double tau_xy = (mu * (d_u_d_y + d_v_d_x));
    const double tau_xz = (mu * (d_u_d_z + d_w_d_x));
    const double tau_yx = tau_xy;
    const double tau_yy = (6.666666666666666e-01 * mu * (-d_u_d_x - d_w_d_z + (2 * d_v_d_y)));
    const double tau_yz = (mu * (d_v_d_z + d_w_d_y));
    const double tau_zx = tau_xz;
    const double tau_zy = tau_yz;
    const double tau_zz = (6.666666666666666e-01 * mu * (-d_u_d_x - d_v_d_y + (2 * d_w_d_z)));
    const double Q_e = ((((d_rho_d_x * u) + (d_u_d_x * rho)) * e_t) + (((d_rho_d_y * v) + (d_v_d_y * rho)) * e_t) + (((d_rho_d_z * w) + (d_w_d_z * rho)) * e_t) - ((d_tau_xx_d_x * u) + (d_u_d_x * tau_xx)) - ((d_tau_xy_d_x * v) + (d_v_d_x * tau_xy)) - ((d_tau_xz_d_x * w) + (d_w_d_x * tau_xz)) - ((d_tau_yx_d_y * u) + (d_u_d_y * tau_yx)) - ((d_tau_yy_d_y * v) + (d_v_d_y * tau_yy)) - ((d_tau_yz_d_y * w) + (d_w_d_y * tau_yz)) - ((d_tau_zx_d_z * u) + (d_u_d_z * tau_zx)) - ((d_tau_zy_d_z * v) + (d_v_d_z * tau_zy)) - ((d_tau_zz_d_z * w) + (d_w_d_z * tau_zz)) + (d_e_t_d_x * rho * u) + (d_e_t_d_y * rho * v) + (d_e_t_d_z * rho * w) + (d_p_d_x * u) + (d_p_d_y * v) + (d_p_d_z * w) + (d_u_d_x * p) + (d_v_d_y * p) + (d_w_d_z * p) + d_q_x_d_x + d_q_y_d_y + d_q_z_d_z);
    */
    // old version
    /*
    const double unique0 = pow(L,(-2));
    const double unique1 = pow(PI,2);
    const double unique4 = pow(L,(-1));
    const double unique102 = (PI * a_rhoz * z * unique4);
    const double unique101 = sin(unique102);
    const double unique105 = (PI * a_px * x * unique4);
    const double unique104 = sin(unique105);
    const double unique107 = (PI * a_py * y * unique4);
    const double unique106 = cos(unique107);
    const double unique111 = (PI * a_pz * z * unique4);
    const double unique110 = sin(unique111);
    const double unique113 = (PI * a_uy * y * unique4);
    const double unique112 = sin(unique113);
    const double unique117 = (PI * a_wy * y * unique4);
    const double unique116 = cos(unique117);
    const double unique3 = (PI * a_ux * x * unique4);
    const double unique122 = sin(unique3);
    const double unique127 = (PI * a_uz * z * unique4);
    const double unique126 = cos(unique127);
    const double unique133 = (PI * a_wx * x * unique4);
    const double unique132 = sin(unique133);
    const double unique82 = (PI * a_vy * y * unique4);
    const double unique152 = sin(unique82);
    const double unique151 = (v_y * unique152);
    const double unique86 = (PI * a_vx * x * unique4);
    const double unique85 = cos(unique86);
    const double unique84 = (v_x * unique85);
    const double unique91 = (PI * a_vz * z * unique4);
    const double unique90 = sin(unique91);
    const double unique89 = (v_z * unique90);
    const double unique150 = (v_0 + unique151 + unique84 + unique89);
    const double unique96 = (PI * a_rhox * x * unique4);
    const double unique153 = cos(unique96);
    const double unique99 = (PI * a_rhoy * y * unique4);
    const double unique154 = sin(unique99);
    const double unique155 = cos(unique102);
    const double unique156 = sin(unique86);
    const double unique157 = cos(unique91);
    const double unique121 = (u_x * unique122);
    const double unique125 = (u_z * unique126);
    const double unique162 = cos(unique113);
    const double unique161 = (u_y * unique162);
    const double unique160 = (u_0 + unique121 + unique125 + unique161);
    const double unique159 = pow(unique160,2);
    const double unique131 = (w_x * unique132);
    const double unique167 = sin(unique117);
    const double unique166 = (w_y * unique167);
    const double unique80 = (PI * a_wz * z * unique4);
    const double unique169 = cos(unique80);
    const double unique168 = (w_z * unique169);
    const double unique165 = (w_0 + unique131 + unique166 + unique168);
    const double unique164 = pow(unique165,2);
    const double unique170 = sin(unique127);
    const double unique174 = cos(unique105);
    const double unique173 = (p_x * unique174);
    const double unique176 = sin(unique107);
    const double unique175 = (p_y * unique176);
    const double unique178 = cos(unique111);
    const double unique177 = (p_z * unique178);
    const double unique172 = (p_0 + unique173 + unique175 + unique177);
    const double unique66 = (-1 + Gamma);
    const double unique65 = pow(unique66,(-1));
    const double unique171 = (Gamma * unique172 * unique4 * unique65);
    const double unique100 = (rho_z * unique101);
    const double unique95 = sin(unique96);
    const double unique94 = (rho_x * unique95);
    const double unique98 = cos(unique99);
    const double unique97 = (rho_y * unique98);
    const double unique93 = (rho_0 + unique100 + unique94 + unique97);
    const double unique179 = pow(unique93,(-3));
    const double unique180 = cos(unique133);
    const double unique181 = pow(unique93,(-1));
    const double unique183 = pow(unique150,2);
    const double unique182 = (unique159 + unique164 + unique183);
    const double unique2 = cos(unique3);
    const double unique29 = pow(R,(-1));
    const double unique79 = sin(unique80);
    const double unique81 = cos(unique82);
    const double unique92 = pow(unique93,(-2));
    const double Q_e = (-(1.333333333333333e+00 * a_ux * a_wz * mu * u_x * w_z * unique0 * unique1 * unique2 * unique79) - (1.333333333333333e+00 * a_vy * a_wz * mu * v_y * w_z * unique0 * unique1 * unique79 * unique81) - (1.333333333333333e+00 * mu * v_y * (-(unique150 * unique152) + (v_y * pow(unique81,2))) * pow(a_vy,2) * unique0 * unique1) - (2 * a_px * a_rhox * k * p_x * rho_x * unique0 * unique1 * unique104 * unique153 * unique29 * unique92) - (2 * a_py * a_rhoy * k * p_y * rho_y * unique0 * unique1 * unique106 * unique154 * unique29 * unique92) - (2 * a_pz * a_rhoz * k * p_z * rho_z * unique0 * unique1 * unique110 * unique155 * unique29 * unique92) - (2 * a_uy * a_vx * mu * u_y * v_x * unique0 * unique1 * unique112 * unique156) - (2 * a_vz * a_wy * mu * v_z * w_y * unique0 * unique1 * unique116 * unique157) - (5.000000000000000e-01 * PI * a_rhoy * rho_y * unique150 * unique154 * unique182 * unique4) - (Gamma * PI * a_px * p_x * unique104 * unique160 * unique4 * unique65) - (Gamma * PI * a_pz * p_z * unique110 * unique165 * unique4 * unique65) - (PI * a_uy * u_y * unique112 * unique150 * unique160 * unique4 * unique93) - (PI * a_uz * u_z * unique160 * unique165 * unique170 * unique4 * unique93) - (a_vx * a_vy * v_x * v_y * ((5.000000000000000e-01 * ((3 * unique183) + unique159 + unique164) * unique4 * unique93) + unique171) * unique1 * unique150 * unique156 * unique160 * unique4 * unique81 * unique93) - (k * rho_x * ((2 * rho_x * pow(unique153,2)) + (unique93 * unique95)) * pow(a_rhox,2) * unique0 * unique1 * unique172 * unique179 * unique29) - (k * rho_y * ((2 * rho_y * pow(unique154,2)) + (unique93 * unique98)) * pow(a_rhoy,2) * unique0 * unique1 * unique172 * unique179 * unique29) - (k * rho_z * ((2 * rho_z * pow(unique155,2)) + (unique101 * unique93)) * pow(a_rhoz,2) * unique0 * unique1 * unique172 * unique179 * unique29) - (mu * v_x * (-(unique150 * unique85) + (v_x * pow(unique156,2))) * pow(a_vx,2) * unique0 * unique1) - (mu * v_z * (-(unique150 * unique90) + (v_z * pow(unique157,2))) * pow(a_vz,2) * unique0 * unique1) + (1.333333333333333e+00 * a_ux * a_vy * mu * u_x * v_y * unique0 * unique1 * unique2 * unique81) + (1.333333333333333e+00 * mu * u_x * (-(u_x * pow(unique2,2)) + (unique122 * unique160)) * pow(a_ux,2) * unique0 * unique1) + (1.333333333333333e+00 * mu * w_z * (-(w_z * pow(unique79,2)) + (unique165 * unique169)) * pow(a_wz,2) * unique0 * unique1) + (2 * a_uz * a_wx * mu * u_z * w_x * unique0 * unique1 * unique170 * unique180) + (5.000000000000000e-01 * PI * a_rhox * rho_x * unique153 * unique160 * unique182 * unique4) + (5.000000000000000e-01 * a_rhoz * a_ux * rho_z * u_x * ((5.000000000000000e-01 * ((3 * unique159) + unique164 + unique183) * unique4 * unique93) + unique171) * unique1 * unique155 * unique165 * unique182 * unique2 * unique4) + (Gamma * PI * a_py * p_y * unique106 * unique150 * unique4 * unique65) + (PI * a_vz * v_z * unique150 * unique157 * unique165 * unique4 * unique93) + (PI * a_wx * w_x * unique160 * unique165 * unique180 * unique4 * unique93) + (a_wy * a_wz * w_y * w_z * (-(5.000000000000000e-01 * ((3 * unique164) + unique159 + unique183) * unique4 * unique93) - unique171) * unique1 * unique116 * unique150 * unique165 * unique4 * unique79 * unique93) + (k * p_x * pow(a_px,2) * unique0 * unique1 * unique174 * unique181 * unique29) + (k * p_y * pow(a_py,2) * unique0 * unique1 * unique176 * unique181 * unique29) + (k * p_z * pow(a_pz,2) * unique0 * unique1 * unique178 * unique181 * unique29) + (mu * u_y * (-(u_y * pow(unique112,2)) + (unique160 * unique162)) * pow(a_uy,2) * unique0 * unique1) + (mu * u_z * (-(u_z * pow(unique170,2)) + (unique126 * unique160)) * pow(a_uz,2) * unique0 * unique1) + (mu * w_x * (-(w_x * pow(unique180,2)) + (unique132 * unique165)) * pow(a_wx,2) * unique0 * unique1) + (mu * w_y * (-(w_y * pow(unique116,2)) + (unique165 * unique167)) * pow(a_wy,2) * unique0 * unique1));
    */
    return(Q_e);
}

double SourceQ_rho (
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
  double Q_rho;
  Q_rho = rho_x * cos(a_rhox * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_rhoy * PI / L + rho_z * cos(a_rhoz * PI * z / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_rhoz * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * a_ux * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * a_vy * PI / L - w_z * sin(a_wz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * a_wz * PI / L;
  return(Q_rho);
}
