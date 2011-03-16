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
// $Id: heat.cpp 18162 2011-03-01 05:23:07Z nick $
//
// heat.cpp: These are the MASA class member functions and constructors
//          For Radiation
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h> 

using namespace MASA;

/* ------------------------------------------------
 *
 *   radiation_integrated_intensity
 *
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
MASA::fans_sa_transient_d_finite<Scalar>::fans_sa_transient_d_finite()
{
    this->mmsname = "fans_sa_transient_d_finite";
    this->dimension=1;

    this->register_var("cb1",&cb1);   
    this->init_var();
  
}//done with constructor

template <typename Scalar>
int MASA::fans_sa_transient_d_finite<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("cb1",1.4);
  return err;

}

template <typename Scalar>
Scalar MASA::fans_sa_transient_d_finite<Scalar>::eval_q_u(Scalar x)
{
  Scalar Q_u;
  Q_u = (cb1 * x);
  return Q_u;
}

// ----------------------------------------
//   Manufactured Analytical Solutions
// ----------------------------------------

template <typename Scalar>
Scalar MASA::fans_sa_transient_d_finite<Scalar>::eval_exact_u(Scalar x)
{
  Scalar exact_u;
  exact_u = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  return exact_u;
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_d_finite<Scalar>::eval_exact_v(Scalar x)
{
  Scalar exact_v;
  exact_v = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);
  return exact_v;
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_d_finite<Scalar>::eval_exact_p(Scalar x)
{
  Scalar exact_p;
  exact_p = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L);
  return exact_p;
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_d_finite<Scalar>::eval_exact_rho(Scalar x)
{
  Scalar exact_rho;
  exact_rho = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  return exact_rho;
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_d_finite<Scalar>::eval_exact_nu(Scalar x)
{
  Scalar exact_nu;
  exact_nu = nu_sa_0 + nu_sa_x * cos(a_nusax * pi * x / L) + nu_sa_y * cos(a_nusay * pi * y / L) + nu_sa_t * cos(a_nusat * pi * t / L);
  return exact_nu;
}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::fans_sa_transient_d_finite);
