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
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <math.h>

typedef double Scalar;

#define PI M_PI
//const Scalar PI = std::acos(-1);

Scalar SourceQ_e ( // 40
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
  Scalar Q_e;
  Q_e = -Gamma * p_x * std::sin(a_px * PI * x / L) * (u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L)) * a_px * PI / L / (Gamma - 0.1e1) + Gamma * p_y * std::cos(a_py * PI * y / L) * (v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L)) * a_py * PI / L / (Gamma - 0.1e1) - Gamma * p_z * std::sin(a_pz * PI * z / L) * (w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L)) * a_pz * PI / L / (Gamma - 0.1e1) + (std::pow(u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L), 0.2e1) + std::pow(v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L), 0.2e1) + std::pow(w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L), 0.2e1)) * a_rhox * PI * rho_x * std::cos(a_rhox * PI * x / L) * (u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L)) / L / 0.2e1 - (std::pow(u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L), 0.2e1) + std::pow(v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L), 0.2e1) + std::pow(w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L), 0.2e1)) * a_rhoy * PI * rho_y * std::sin(a_rhoy * PI * y / L) * (v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L)) / L / 0.2e1 + (std::pow(u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L), 0.2e1) + std::pow(v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L), 0.2e1) + std::pow(w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L), 0.2e1)) * a_rhoz * PI * rho_z * std::cos(a_rhoz * PI * z / L) * (w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L)) / L / 0.2e1 + a_ux * PI * u_x * std::cos(a_ux * PI * x / L) * ((0.3e1 * rho_x * std::sin(a_rhox * PI * x / L) + 0.3e1 * rho_y * std::cos(a_rhoy * PI * y / L) + 0.3e1 * rho_z * std::sin(a_rhoz * PI * z / L) + 0.3e1 * rho_0) * std::pow(u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L), 0.2e1) + (rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L) + rho_0) * std::pow(v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L), 0.2e1) + (rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L) + rho_0) * std::pow(w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L), 0.2e1)) / L / 0.2e1 + Gamma * a_ux * PI * u_x * std::cos(a_ux * PI * x / L) * (p_x * std::cos(a_px * PI * x / L) + p_y * std::sin(a_py * PI * y / L) + p_z * std::cos(a_pz * PI * z / L) + p_0) / L / (Gamma - 0.1e1) - (v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L)) * (rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L) + rho_0) * (u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L)) * PI * a_uy * u_y * std::sin(a_uy * PI * y / L) / L - (w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L)) * (rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L) + rho_0) * (u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L)) * PI * a_uz * u_z * std::sin(a_uz * PI * z / L) / L - (v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L)) * (rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L) + rho_0) * (u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L)) * PI * a_vx * v_x * std::sin(a_vx * PI * x / L) / L + a_vy * PI * v_y * std::cos(a_vy * PI * y / L) * ((rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L) + rho_0) * std::pow(u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L), 0.2e1) + (0.3e1 * rho_x * std::sin(a_rhox * PI * x / L) + 0.3e1 * rho_y * std::cos(a_rhoy * PI * y / L) + 0.3e1 * rho_z * std::sin(a_rhoz * PI * z / L) + 0.3e1 * rho_0) * std::pow(v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L), 0.2e1) + (rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L) + rho_0) * std::pow(w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L), 0.2e1)) / L / 0.2e1 + Gamma * a_vy * PI * v_y * std::cos(a_vy * PI * y / L) * (p_x * std::cos(a_px * PI * x / L) + p_y * std::sin(a_py * PI * y / L) + p_z * std::cos(a_pz * PI * z / L) + p_0) / L / (Gamma - 0.1e1) + (w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L)) * (rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L) + rho_0) * (v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L)) * PI * a_vz * v_z * std::cos(a_vz * PI * z / L) / L + (w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L)) * (rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L) + rho_0) * (u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L)) * PI * a_wx * w_x * std::cos(a_wx * PI * x / L) / L + (w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L)) * (rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L) + rho_0) * (v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L)) * PI * a_wy * w_y * std::cos(a_wy * PI * y / L) / L - a_wz * PI * w_z * std::sin(a_wz * PI * z / L) * ((rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L) + rho_0) * std::pow(u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L), 0.2e1) + (rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L) + rho_0) * std::pow(v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L), 0.2e1) + (0.3e1 * rho_x * std::sin(a_rhox * PI * x / L) + 0.3e1 * rho_y * std::cos(a_rhoy * PI * y / L) + 0.3e1 * rho_z * std::sin(a_rhoz * PI * z / L) + 0.3e1 * rho_0) * std::pow(w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L), 0.2e1)) / L / 0.2e1 - Gamma * a_wz * PI * w_z * std::sin(a_wz * PI * z / L) * (p_x * std::cos(a_px * PI * x / L) + p_y * std::sin(a_py * PI * y / L) + p_z * std::cos(a_pz * PI * z / L) + p_0) / L / (Gamma - 0.1e1);
  return(Q_e);
}

Scalar SourceQ_u ( // 38 variables
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
  Scalar Q_u;
  Q_u = -p_x * std::sin(a_px * PI * x / L) * a_px * PI / L + rho_x * std::cos(a_rhox * PI * x / L) * std::pow(u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L), 0.2e1) * a_rhox * PI / L - rho_y * std::sin(a_rhoy * PI * y / L) * (v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L)) * (u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L)) * a_rhoy * PI / L + rho_z * std::cos(a_rhoz * PI * z / L) * (w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L)) * (u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L)) * a_rhoz * PI / L + 0.2e1 * u_x * std::cos(a_ux * PI * x / L) * (rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L)) * (u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L)) * a_ux * PI / L - u_y * std::sin(a_uy * PI * y / L) * (rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L)) * (v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L)) * a_uy * PI / L - u_z * std::sin(a_uz * PI * z / L) * (rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L)) * (w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L)) * a_uz * PI / L + v_y * std::cos(a_vy * PI * y / L) * (rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L)) * (u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L)) * a_vy * PI / L - w_z * std::sin(a_wz * PI * z / L) * (rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L)) * (u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L)) * a_wz * PI / L;
  return(Q_u);
}

Scalar SourceQ_v (
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
  Scalar Q_v;
  Q_v = p_y * std::cos(a_py * PI * y / L) * a_py * PI / L + rho_x * std::cos(a_rhox * PI * x / L) * (v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L)) * (u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L)) * a_rhox * PI / L - rho_y * std::sin(a_rhoy * PI * y / L) * std::pow(v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L), 0.2e1) * a_rhoy * PI / L + rho_z * std::cos(a_rhoz * PI * z / L) * (w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L)) * (v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L)) * a_rhoz * PI / L + u_x * std::cos(a_ux * PI * x / L) * (rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L)) * (v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L)) * a_ux * PI / L - v_x * std::sin(a_vx * PI * x / L) * (rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L)) * (u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L)) * a_vx * PI / L + 0.2e1 * v_y * std::cos(a_vy * PI * y / L) * (rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L)) * (v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L)) * a_vy * PI / L + v_z * std::cos(a_vz * PI * z / L) * (rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L)) * (w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L)) * a_vz * PI / L - w_z * std::sin(a_wz * PI * z / L) * (rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L)) * (v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L)) * a_wz * PI / L;
  return(Q_v);
}

Scalar SourceQ_w (
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
  Scalar Q_w;
  Q_w = -p_z * std::sin(a_pz * PI * z / L) * a_pz * PI / L + rho_x * std::cos(a_rhox * PI * x / L) * (w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L)) * (u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L)) * a_rhox * PI / L - rho_y * std::sin(a_rhoy * PI * y / L) * (w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L)) * (v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L)) * a_rhoy * PI / L + rho_z * std::cos(a_rhoz * PI * z / L) * std::pow(w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L), 0.2e1) * a_rhoz * PI / L + u_x * std::cos(a_ux * PI * x / L) * (rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L)) * (w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L)) * a_ux * PI / L + v_y * std::cos(a_vy * PI * y / L) * (rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L)) * (w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L)) * a_vy * PI / L + w_x * std::cos(a_wx * PI * x / L) * (rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L)) * (u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L)) * a_wx * PI / L + w_y * std::cos(a_wy * PI * y / L) * (rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L)) * (v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L)) * a_wy * PI / L - 0.2e1 * w_z * std::sin(a_wz * PI * z / L) * (rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L)) * (w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L)) * a_wz * PI / L;
  return(Q_w);
}

Scalar SourceQ_rho(
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
  Scalar Q_rho;
  Q_rho = rho_x * std::cos(a_rhox * PI * x / L) * (u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L) + u_z * std::cos(a_uz * PI * z / L)) * a_rhox * PI / L - rho_y * std::sin(a_rhoy * PI * y / L) * (v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L) + v_z * std::sin(a_vz * PI * z / L)) * a_rhoy * PI / L + rho_z * std::cos(a_rhoz * PI * z / L) * (w_0 + w_x * std::sin(a_wx * PI * x / L) + w_y * std::sin(a_wy * PI * y / L) + w_z * std::cos(a_wz * PI * z / L)) * a_rhoz * PI / L + u_x * std::cos(a_ux * PI * x / L) * (rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L)) * a_ux * PI / L + v_y * std::cos(a_vy * PI * y / L) * (rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L)) * a_vy * PI / L - w_z * std::sin(a_wz * PI * z / L) * (rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L) + rho_z * std::sin(a_rhoz * PI * z / L)) * a_wz * PI / L;
  return(Q_rho);
}
