// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012,2013 The PECOS Development Team
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
// axi_euler.cpp: program that tests axisymmetric euler
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <tests.h>
#include <cmath>

using namespace MASA;
using namespace std;

template<typename Scalar>
Scalar anQ_p(Scalar r,Scalar z,Scalar p_0,Scalar p_r,Scalar p_z,Scalar rho_0,Scalar rho_r,Scalar u_r,Scalar w_0,Scalar w_r,Scalar a_pr,Scalar a_pz,Scalar a_rhor,Scalar a_rhoz,Scalar a_ur,Scalar a_uz,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L,Scalar Gamma)
{
  Scalar exact_p = p_0 + p_r * std::sin(a_pr * pi * r / L) + p_z * std::cos(a_pz * pi * z / L);
  return exact_p;
}

template<typename Scalar>
Scalar anQ_u (Scalar r,Scalar z,Scalar p_0,Scalar p_r,Scalar rho_0,Scalar rho_r,Scalar u_r,Scalar u_z,Scalar w_0,Scalar w_r,Scalar a_pr,Scalar a_pz,Scalar a_rhor,Scalar a_rhoz,Scalar a_ur,Scalar a_uz,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L,Scalar Gamma)
{
  Scalar exact_u = u_r * u_z * (std::cos(a_ur * pi * r / L) - Scalar(0.1e1)) * std::sin(a_uz * pi * z / L);
  return exact_u;
}

template<typename Scalar>
Scalar anQ_w (Scalar r,Scalar z,Scalar w_0,Scalar w_r,Scalar w_z,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L)
{
  Scalar exact_w = w_0 + w_r * std::cos(a_wr * pi * r / L) + w_z * std::sin(a_wz * pi * z / L);
  return exact_w;
}

template<typename Scalar>
Scalar anQ_rho (Scalar r,Scalar z,Scalar p_0,Scalar p_r,Scalar rho_0,Scalar rho_r,Scalar rho_z,Scalar u_r,Scalar w_0,Scalar w_r,Scalar a_pr,Scalar a_pz,Scalar a_rhor,Scalar a_rhoz,Scalar a_ur,Scalar a_uz,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L,Scalar Gamma)
{
  Scalar exact_rho = rho_0 + rho_r * std::cos(a_rhor * pi * r / L) + rho_z * std::sin(a_rhoz * pi * z / L);
  return exact_rho;
}

// ----------------------------------------
//   Source Terms
// ----------------------------------------

template<typename Scalar>
Scalar SourceQ_e(Scalar r,Scalar z,Scalar p_0,Scalar p_r,Scalar p_z,Scalar rho_0,Scalar rho_r,Scalar rho_z,Scalar u_r,Scalar u_z,Scalar w_0,Scalar w_r,Scalar w_z,Scalar a_pr,Scalar a_pz,Scalar a_rhor,Scalar a_rhoz,Scalar a_ur,Scalar a_uz,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L,Scalar Gamma)
{
  Scalar Q_e;
  Scalar RHO;
  Scalar P;
  Scalar U;
  Scalar W;
  RHO = rho_0 + rho_r * std::cos(a_rhor * pi * r / L) + rho_z * std::sin(a_rhoz * pi * z / L);
  P = p_0 + p_r * std::sin(a_pr * pi * r / L) + p_z * std::cos(a_pz * pi * z / L);
  U = u_r * u_z * (std::cos(a_ur * pi * r / L) - 0.1e1) * std::sin(a_uz * pi * z / L);
  W = w_0 + w_r * std::cos(a_wr * pi * r / L) + w_z * std::sin(a_wz * pi * z / L);
  Q_e = -Gamma * a_ur * pi * u_r * u_z * P * std::sin(a_ur * pi * r / L) * std::sin(a_uz * pi * z / L) / (Gamma - 0.1e1) / L + (std::cos(a_ur * pi * r / L) - 0.1e1) * a_uz * pi * u_r * u_z * RHO * U * W * std::cos(a_uz * pi * z / L) / L - (0.3e1 * U * U + W * W) * a_ur * pi * u_r * u_z * RHO * std::sin(a_ur * pi * r / L) * std::sin(a_uz * pi * z / L) / L / 0.2e1 - a_wr * pi * w_r * RHO * U * W * std::sin(a_wr * pi * r / L) / L + Gamma * a_pr * pi * p_r * U * std::cos(a_pr * pi * r / L) / (Gamma - 0.1e1) / L - Gamma * a_pz * pi * p_z * W * std::sin(a_pz * pi * z / L) / (Gamma - 0.1e1) / L + Gamma * a_wz * pi * w_z * P * std::cos(a_wz * pi * z / L) / (Gamma - 0.1e1) / L - (U * U + W * W) * a_rhor * pi * rho_r * U * std::sin(a_rhor * pi * r / L) / L / 0.2e1 + (U * U + W * W) * a_rhoz * pi * rho_z * W * std::cos(a_rhoz * pi * z / L) / L / 0.2e1 + (U * U + 0.3e1 * W * W) * a_wz * pi * w_z * RHO * std::cos(a_wz * pi * z / L) / L / 0.2e1 + Gamma * P * U / (Gamma - 0.1e1) / r + (U * U + W * W) * RHO * U / r / 0.2e1;
  return(Q_e);
}

template<typename Scalar>
Scalar SourceQ_u(Scalar r,Scalar z,Scalar p_0,Scalar p_r,Scalar rho_0,Scalar rho_r,Scalar rho_z,Scalar u_r,Scalar u_z,Scalar w_0,Scalar w_r,Scalar w_z,Scalar a_pr,Scalar a_pz,Scalar a_rhor,Scalar a_rhoz,Scalar a_ur,Scalar a_uz,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L,Scalar Gamma)
{
  Scalar Q_u;
  Scalar RHO;
  Scalar U;
  Scalar W;
  RHO = rho_0 + rho_r * std::cos(a_rhor * pi * r / L) + rho_z * std::sin(a_rhoz * pi * z / L);
  U = u_r * u_z * (std::cos(a_ur * pi * r / L) - 0.1e1) * std::sin(a_uz * pi * z / L);
  W = w_0 + w_r * std::cos(a_wr * pi * r / L) + w_z * std::sin(a_wz * pi * z / L);
  Q_u = (std::cos(a_ur * pi * r / L) - 0.1e1) * a_uz * pi * u_r * u_z * RHO * W * std::cos(a_uz * pi * z / L) / L - a_rhor * pi * rho_r * U * U * std::sin(a_rhor * pi * r / L) / L + a_rhoz * pi * rho_z * U * W * std::cos(a_rhoz * pi * z / L) / L + a_pr * pi * p_r * std::cos(a_pr * pi * r / L) / L - (0.2e1 * a_ur * u_r * u_z * std::sin(a_ur * pi * r / L) * std::sin(a_uz * pi * z / L) - a_wz * w_z * std::cos(a_wz * pi * z / L)) * pi * RHO * U / L + RHO * U * U / r;
  return(Q_u);
}

template<typename Scalar>
Scalar SourceQ_w(Scalar r,Scalar z,Scalar p_0,Scalar p_r,Scalar p_z,Scalar rho_0,Scalar rho_r,Scalar rho_z,Scalar u_r,Scalar u_z,Scalar w_0,Scalar w_r,Scalar w_z,Scalar a_pr,Scalar a_pz,Scalar a_rhor,Scalar a_rhoz,Scalar a_ur,Scalar a_uz,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L,Scalar Gamma)
{
  Scalar Q_w;
  Scalar RHO;
  Scalar U;
  Scalar W;
  RHO = rho_0 + rho_r * std::cos(a_rhor * pi * r / L) + rho_z * std::sin(a_rhoz * pi * z / L);
  U = u_r * u_z * (std::cos(a_ur * pi * r / L) - 0.1e1) * std::sin(a_uz * pi * z / L);
  W = w_0 + w_r * std::cos(a_wr * pi * r / L) + w_z * std::sin(a_wz * pi * z / L);
  Q_w = -a_rhor * pi * rho_r * U * W * std::sin(a_rhor * pi * r / L) / L + a_rhoz * pi * rho_z * W * W * std::cos(a_rhoz * pi * z / L) / L - a_wr * pi * w_r * RHO * U * std::sin(a_wr * pi * r / L) / L - a_pz * pi * p_z * std::sin(a_pz * pi * z / L) / L - (a_ur * u_r * u_z * std::sin(a_ur * pi * r / L) * std::sin(a_uz * pi * z / L) - 0.2e1 * a_wz * w_z * std::cos(a_wz * pi * z / L)) * pi * RHO * W / L + RHO * U * W / r;
  return(Q_w);
}

template<typename Scalar>
Scalar SourceQ_rho(Scalar r,Scalar z,Scalar p_0,Scalar p_r,Scalar rho_0,Scalar rho_r,Scalar rho_z,Scalar u_r,Scalar u_z,Scalar w_0,Scalar w_r,Scalar w_z,Scalar a_pr,Scalar a_pz,Scalar a_rhor,Scalar a_rhoz,Scalar a_ur,Scalar a_uz,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L,Scalar Gamma)
{
  Scalar Q_rho;
  Scalar RHO;
  Scalar U;
  Scalar W;
  RHO = rho_0 + rho_r * std::cos(a_rhor * pi * r / L) + rho_z * std::sin(a_rhoz * pi * z / L);
  U = u_r * u_z * (std::cos(a_ur * pi * r / L) - 0.1e1) * std::sin(a_uz * pi * z / L);
  W = w_0 + w_r * std::cos(a_wr * pi * r / L) + w_z * std::sin(a_wz * pi * z / L);
  Q_rho = -a_rhor * pi * rho_r * U * std::sin(a_rhor * pi * r / L) / L + a_rhoz * pi * rho_z * W * std::cos(a_rhoz * pi * z / L) / L - (a_ur * u_r * u_z * std::sin(a_ur * pi * r / L) * std::sin(a_uz * pi * z / L) - a_wz * w_z * std::cos(a_wz * pi * z / L)) * pi * RHO / L + RHO * U / r;
  return(Q_rho);
}

template<typename Scalar>
int run_regression()
{
  //variables
  Scalar p_0;
  Scalar p_r;
  Scalar p_z;
  Scalar rho_0;
  Scalar rho_r;
  Scalar rho_z;
  Scalar u_r;
  Scalar u_z;
  Scalar w_0;
  Scalar w_r;
  Scalar w_z;
  Scalar a_pr;
  Scalar a_pz;
  Scalar a_rhor;
  Scalar a_rhoz;
  Scalar a_ur;
  Scalar a_uz;
  Scalar a_wr;
  Scalar a_wz;
  Scalar L;
  Scalar Gamma;

  // parameters
  Scalar r;
  Scalar z;

  // solutions
  Scalar ufield,ufield2,ufield3;
  Scalar wfield,wfield2,wfield3;
  Scalar efield,efield2,efield3;
  Scalar rho,rho2,rho3;

  Scalar exact_u,exact_u2,exact_u3;
  Scalar exact_w,exact_w2,exact_w3;
  Scalar exact_p,exact_p2,exact_p3;
  Scalar exact_rho,exact_rho2,exact_rho3;

  const Scalar pi = std::acos(Scalar(-1));

  // initalize
  int nx = 115;  // number of points
  int ny = 68;
  int lx=3;     // length
  int ly=1;

  Scalar dx=Scalar(lx)/Scalar(nx);
  Scalar dy=Scalar(ly)/Scalar(ny);

  masa_init<Scalar>("axisymmetric_euler","axisymmetric_euler");

  // set params
  masa_init_param<Scalar>();

  // get vars
  p_0    = masa_get_param<Scalar>("p_0");
  p_r    = masa_get_param<Scalar>("p_r");
  p_z    = masa_get_param<Scalar>("p_z");
  rho_0  = masa_get_param<Scalar>("rho_0");
  rho_r  = masa_get_param<Scalar>("rho_r");
  rho_z  = masa_get_param<Scalar>("rho_z");
  u_r    = masa_get_param<Scalar>("u_r");
  u_z    = masa_get_param<Scalar>("u_z");
  w_0    = masa_get_param<Scalar>("w_0");
  w_r    = masa_get_param<Scalar>("w_r");
  w_z    = masa_get_param<Scalar>("w_z");
  a_pr   = masa_get_param<Scalar>("a_pr");
  a_pz   = masa_get_param<Scalar>("a_pz");
  a_rhor = masa_get_param<Scalar>("a_rhor");
  a_rhoz = masa_get_param<Scalar>("a_rhoz");
  a_ur   = masa_get_param<Scalar>("a_ur");
  a_uz   = masa_get_param<Scalar>("a_uz");
  a_wr   = masa_get_param<Scalar>("a_wr");
  a_wz   = masa_get_param<Scalar>("a_wz");
  L      = masa_get_param<Scalar>("L");
  Gamma  = masa_get_param<Scalar>("Gamma");

  // check that all terms have been initialized
  int err = masa_sanity_check<Scalar>();
  if(err != 0)
    {
      cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }

  // evaluate source terms (2D)
  for(int i=1;i<nx;i++)    // this is the radial term -- thus, do not start at 0!
    for(int j=1;j<ny;j++)  // z component
      {
	r=i*dx;
	z=j*dy;

	//evalulate source terms
	ufield = masa_eval_source_rho_u<Scalar>  (r,z);
	wfield = masa_eval_source_rho_w<Scalar>  (r,z);
	efield = masa_eval_source_rho_e<Scalar>  (r,z);
	rho    = masa_eval_source_rho<Scalar>(r,z);

	//evaluate analytical terms
	exact_u = masa_eval_exact_u<Scalar>        (r,z);
	exact_w = masa_eval_exact_w<Scalar>        (r,z);
	exact_p = masa_eval_exact_p<Scalar>        (r,z);
	exact_rho = masa_eval_exact_rho<Scalar>    (r,z);

	// check against maple
	ufield2 = SourceQ_u   (r, z, p_0, p_r, rho_0, rho_r, rho_z, u_r, u_z, w_0, w_r, w_z, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);
	wfield2 = SourceQ_w   (r, z, p_0, p_r, p_z, rho_0, rho_r, rho_z, u_r, u_z, w_0, w_r, w_z, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);
	rho2    = SourceQ_rho (r, z, p_0, p_r, rho_0, rho_r, rho_z, u_r, u_z, w_0, w_r, w_z, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);
	efield2 = SourceQ_e   (r, z, p_0, p_r, p_z, rho_0, rho_r, rho_z, u_r, u_z, w_0, w_r, w_z, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);

	exact_u2   = anQ_u   (r, z, p_0, p_r, rho_0, rho_r, u_r, u_z, w_0, w_r, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);
	exact_rho2 = anQ_rho (r, z, p_0, p_r, rho_0, rho_r, rho_z, u_r, w_0, w_r, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);
	exact_p2   = anQ_p   (r, z, p_0, p_r, p_z, rho_0, rho_r, u_r, w_0, w_r, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);
	exact_w2   = anQ_w   (r, z, w_0, w_r, w_z, a_wr, a_wz, pi, L);

	// test the result is roughly zero
	// choose between abs and rel error
#ifdef MASA_STRICT_REGRESSION

	ufield3 = std::abs(ufield-ufield2);
	wfield3 = std::abs(wfield-wfield2);
	efield3 = std::abs(efield-efield2);
	rho3    = std::abs(rho-rho2);

	exact_u3   = std::abs(exact_u-exact_u2);
	exact_w3   = std::abs(exact_w-exact_w2);
	exact_rho3 = std::abs(exact_rho-exact_rho2);
	exact_p3   = std::abs(exact_p-exact_p2);

#else

	ufield3 = std::abs(ufield-ufield2)/std::abs(ufield2);
	wfield3 = std::abs(wfield-wfield2)/std::abs(wfield2);
	efield3 = std::abs(efield-efield2)/std::abs(efield2);
	rho3    = std::abs(rho-rho2)/std::abs(rho2);

	exact_u3   = std::abs(exact_u-exact_u2)/std::abs(exact_u2);
	exact_w3   = std::abs(exact_w-exact_w2)/std::abs(exact_w2);
	exact_rho3 = std::abs(exact_rho-exact_rho2)/std::abs(exact_rho2);
	exact_p3   = std::abs(exact_p-exact_p2)/std::abs(exact_p2);

#endif

	threshcheck(ufield3);
	threshcheck(wfield3);
	threshcheck(efield3);
	threshcheck(rho3);

	threshcheck(exact_u3);
	threshcheck(exact_w3);
	threshcheck(exact_rho3);
	threshcheck(exact_p3);

      } // done iterating

  // tests passed
  return 0;
}

int main()
{
  int err=0;

  err += run_regression<double>();
  err += run_regression<long double>();

  return err;
}
