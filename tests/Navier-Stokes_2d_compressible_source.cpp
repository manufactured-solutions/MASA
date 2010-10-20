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

typedef double Scalar;

const Scalar PI = acos(-1);

Scalar SourceQ_u (
  Scalar x,
  Scalar y,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar a_px,
  Scalar a_py,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_vx,
  Scalar a_vy,
  Scalar mu,
  Scalar L,
  Scalar R,
  Scalar k)
{
  Scalar Q_u;
  Q_u = 0.4e1 / 0.3e1 * mu * u_x * sin(a_ux * PI * x / L) * a_ux * a_ux * PI * PI * pow(L, -0.2e1) + mu * u_y * cos(a_uy * PI * y / L) * a_uy * a_uy * PI * PI * pow(L, -0.2e1) - p_x * sin(a_px * PI * x / L) * a_px * PI / L + rho_x * cos(a_rhox * PI * x / L) * pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L), 0.2e1) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_rhoy * PI / L + 0.2e1 * u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_ux * PI / L - u_y * sin(a_uy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_uy * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_vy * PI / L;
  return(Q_u);
}

Scalar SourceQ_v (
  Scalar x,
  Scalar y,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar a_px,
  Scalar a_py,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_vx,
  Scalar a_vy,
  Scalar mu,
  Scalar L,
  Scalar R,
  Scalar k)
{
  Scalar Q_v;
  Q_v = mu * v_x * cos(a_vx * PI * x / L) * a_vx * a_vx * PI * PI * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * mu * v_y * sin(a_vy * PI * y / L) * a_vy * a_vy * PI * PI * pow(L, -0.2e1) + p_y * cos(a_py * PI * y / L) * a_py * PI / L + rho_x * cos(a_rhox * PI * x / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L), 0.2e1) * a_rhoy * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_ux * PI / L - v_x * sin(a_vx * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_vx * PI / L + 0.2e1 * v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_vy * PI / L;
  return(Q_v);
}

#include <math.h>

Scalar SourceQ_e (
  Scalar x,
  Scalar y,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar a_px,
  Scalar a_py,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_vx,
  Scalar a_vy,
  Scalar Gamma,
  Scalar mu,
  Scalar L,
  Scalar R,
  Scalar k)
{
  Scalar Q_e;
  Q_e = -(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * (pow(u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1) + pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0, 0.2e1)) * rho_y * sin(a_rhoy * PI * y / L) * a_rhoy * PI / L / 0.2e1 + (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * (pow(u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1) + pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0, 0.2e1)) * rho_x * cos(a_rhox * PI * x / L) * a_rhox * PI / L / 0.2e1 + 0.4e1 / 0.3e1 * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * mu * v_y * sin(a_vy * PI * y / L) * a_vy * a_vy * PI * PI * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * mu * v_y * v_y * pow(cos(a_vy * PI * y / L), 0.2e1) * a_vy * a_vy * PI * PI * pow(L, -0.2e1) - mu * v_x * v_x * pow(sin(a_vx * PI * x / L), 0.2e1) * a_vx * a_vx * PI * PI * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * mu * u_x * u_x * pow(cos(a_ux * PI * x / L), 0.2e1) * a_ux * a_ux * PI * PI * pow(L, -0.2e1) - mu * u_y * u_y * pow(sin(a_uy * PI * y / L), 0.2e1) * a_uy * a_uy * PI * PI * pow(L, -0.2e1) + (Gamma * (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_0) / (Gamma - 0.1e1) + (pow(u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1) / 0.2e1 + 0.3e1 / 0.2e1 * pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0, 0.2e1)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0)) * v_y * cos(a_vy * PI * y / L) * a_vy * PI / L + (Gamma * (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_0) / (Gamma - 0.1e1) + (0.3e1 / 0.2e1 * pow(u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1) + pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0, 0.2e1) / 0.2e1) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0)) * u_x * cos(a_ux * PI * x / L) * a_ux * PI / L + (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * mu * v_x * cos(a_vx * PI * x / L) * a_vx * a_vx * PI * PI * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * mu * u_x * sin(a_ux * PI * x / L) * a_ux * a_ux * PI * PI * pow(L, -0.2e1) + (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * mu * u_y * cos(a_uy * PI * y / L) * a_uy * a_uy * PI * PI * pow(L, -0.2e1) - (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * u_y * sin(a_uy * PI * y / L) * a_uy * PI / L - (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_0) * rho_x * k * sin(a_rhox * PI * x / L) * a_rhox * a_rhox * PI * PI * pow(rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0, -0.2e1) * pow(L, -0.2e1) / R - (0.2e1 * p_x * cos(a_px * PI * x / L) + 0.2e1 * p_y * sin(a_py * PI * y / L) + 0.2e1 * p_0) * rho_x * rho_x * k * pow(cos(a_rhox * PI * x / L), 0.2e1) * a_rhox * a_rhox * PI * PI * pow(rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0, -0.3e1) * pow(L, -0.2e1) / R - (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_0) * rho_y * k * cos(a_rhoy * PI * y / L) * a_rhoy * a_rhoy * PI * PI * pow(rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0, -0.2e1) * pow(L, -0.2e1) / R - (0.2e1 * p_x * cos(a_px * PI * x / L) + 0.2e1 * p_y * sin(a_py * PI * y / L) + 0.2e1 * p_0) * rho_y * rho_y * k * pow(sin(a_rhoy * PI * y / L), 0.2e1) * a_rhoy * a_rhoy * PI * PI * pow(rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0, -0.3e1) * pow(L, -0.2e1) / R + 0.4e1 / 0.3e1 * mu * u_x * v_y * cos(a_ux * PI * x / L) * cos(a_vy * PI * y / L) * a_ux * a_vy * PI * PI * pow(L, -0.2e1) - 0.2e1 * mu * u_y * v_x * sin(a_uy * PI * y / L) * sin(a_vx * PI * x / L) * a_uy * a_vx * PI * PI * pow(L, -0.2e1) - 0.2e1 * k * p_x * rho_x * cos(a_rhox * PI * x / L) * sin(a_px * PI * x / L) * a_px * a_rhox * PI * PI * pow(rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0, -0.2e1) * pow(L, -0.2e1) / R - 0.2e1 * k * p_y * rho_y * cos(a_py * PI * y / L) * sin(a_rhoy * PI * y / L) * a_py * a_rhoy * PI * PI * pow(rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0, -0.2e1) * pow(L, -0.2e1) / R - (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * v_x * sin(a_vx * PI * x / L) * a_vx * PI / L - Gamma * (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * p_x * sin(a_px * PI * x / L) * a_px * PI / (Gamma - 0.1e1) / L + Gamma * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * p_y * cos(a_py * PI * y / L) * a_py * PI / (Gamma - 0.1e1) / L + k * p_x * cos(a_px * PI * x / L) * a_px * a_px * PI * PI / (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * pow(L, -0.2e1) / R + k * p_y * sin(a_py * PI * y / L) * a_py * a_py * PI * PI / (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * pow(L, -0.2e1) / R;
  return(Q_e);
}

Scalar SourceQ_rho (
  Scalar x,
  Scalar y,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar a_px,
  Scalar a_py,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_vx,
  Scalar a_vy,
  Scalar mu,
  Scalar L,
  Scalar R,
  Scalar k)

{
  Scalar Q_rho;
  Q_rho = (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * a_rhox * PI * rho_x * cos(a_rhox * PI * x / L) / L - (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * a_rhoy * PI * rho_y * sin(a_rhoy * PI * y / L) / L + (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * a_ux * PI * u_x * cos(a_ux * PI * x / L) / L + (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * a_vy * PI * v_y * cos(a_vy * PI * y / L) / L;
  return(Q_rho);
}
