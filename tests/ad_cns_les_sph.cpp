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
//
// ad_cns.cpp: test the incompressible navier stokes equations are roughly identical to
//              the source terms generated from AD
//
// $Id: $
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <iostream>
#include <config.h>

#ifdef HAVE_METAPHYSICL

#include <tests.h>
#include "ad_masa.h"

//typedef double RawScalar;
typedef ShadowNumber<double, long double> RawScalar;

const unsigned int NDIM = 3;

typedef DualNumber<RawScalar, NumberVector<NDIM, RawScalar> > FirstDerivType;
typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;

typedef SecondDerivType ADType;

template <std::size_t NDIM, typename Scalar>
double evaluate_q (const NumberVector<NDIM, Scalar>& xyz, const int);

using namespace MASA;

int main(void)
{
  int err = 0;
  int N   = 10; // mesh pts. in x, y, z
  double h = 1.0/N;
  double su, sv, sw, s2u, s2v, s2w, sp, se, s2e, s2p;
  double pnorm, unorm, vnorm, wnorm, enorm;
  double pnorm_max, unorm_max, vnorm_max, wnorm_max, enorm_max;
  double prnorm_max = 1., urnorm_max = 0., vrnorm_max = 0., wrnorm_max, ernorm_max = 0.;

  unorm_max = 0;
  vnorm_max = 0;
  wnorm_max = 0;
  pnorm_max = 0;
  enorm_max = 0;

  const RawScalar xvecinit[] = {1., 0., 0.};
  const RawScalar yvecinit[] = {0., 1., 0.};
  const RawScalar zvecinit[] = {0., 0., 1.};

  const NumberVector<NDIM, RawScalar> xvec(xvecinit);
  const NumberVector<NDIM, RawScalar> yvec(yvecinit);
  const NumberVector<NDIM, RawScalar> zvec(zvecinit);

  // initialize the problem in MASA
  err += masa_init("cnd_les_sph","ad_cns_3d_les_sph");

  // If you are wondering... these are CGS units, roughly
  // corresponding to a Re=1, Mach=1, and Prandtl=1 flow.
  masa_set_param("L", 2.);
  masa_set_param("R", 2870024.971535356);
  masa_set_param("k", 820641949.8052748);
  masa_set_param("Gamma", 1.4);
  masa_set_param("mu", 81.69584963240251);
  masa_set_param("mu_bulk", 81.69584963240251);
  masa_set_param("rho_0", 0.001176528671407439);
  masa_set_param("rho_x", 0.0001176528671407439);
  masa_set_param("rho_y", 1.);
  masa_set_param("rho_z", 1.);
  masa_set_param("u_0", 34719.02199147968);
  masa_set_param("u_x", 3471.902199147968);
  masa_set_param("u_y", 1.);
  masa_set_param("u_z", 1.);
  masa_set_param("v_0", 34719.02199147968);
  masa_set_param("v_x", 1.);
  masa_set_param("v_y", 3471.902199147968);
  masa_set_param("v_z", 1.);
  masa_set_param("w_0", 34719.02199147968);
  masa_set_param("w_x", 1.);
  masa_set_param("w_y", 1.);
  masa_set_param("w_z", 3471.902199147968);
  masa_set_param("p_0", 1013000.);
  masa_set_param("p_x", 202600.);
  masa_set_param("p_y", 1.);
  masa_set_param("p_z", 1.);
  masa_set_param("a_rhox", 2.1);
  masa_set_param("a_rhoy", 4.2);
  masa_set_param("a_rhoz", 6.0);
  masa_set_param("a_ux", 2.3);
  masa_set_param("a_uy", 4.4);
  masa_set_param("a_uz", 6.1);
  masa_set_param("a_vx", 4.7);
  masa_set_param("a_vy", 2.2);
  masa_set_param("a_vz", 6.6);
  masa_set_param("a_wx", 6.2);
  masa_set_param("a_wy", 4.8);
  masa_set_param("a_wz", 2.1);
  masa_set_param("a_px", 6.9);
  masa_set_param("a_py", 2.3);
  masa_set_param("a_pz", 4.4);
  masa_set_param("Cs", 0.16);
  masa_set_param("CI", 0.09);
  masa_set_param("PrT", 0.7);
  masa_set_param("deltabar", h);

  // call the sanity check routine
  err += masa_sanity_check();

  // we first set up the DualNumbers that correspond to independent
  // variables, spatial coordinates x, y, z.
  NumberVector<NDIM, ADType> xyz;
  xyz[0] = ADType(1., xvec);
  xyz[1] = ADType(1., yvec);
  xyz[2] = ADType(1., zvec);

  double x = 0, y = 0, z = 0;
  for (int i=0; i != N+1; ++i){
    x = i*h;
    xyz[0] = ADType(x, xvec);

    for (int j=0; j != N+1; ++j){
      y = j*h;
      xyz[1] = ADType(y, yvec);

      for (int k=0; k != N+1; ++k){
	z = k*h;
	xyz[2] = ADType(k*h, zvec);

	// evaluate masa source terms
	su = masa_eval_source_rho_u<double>(x,y,z);
	sv = masa_eval_source_rho_v<double>(x,y,z);
	sw = masa_eval_source_rho_w<double>(x,y,z);
	sp = masa_eval_source_rho  <double>(x,y,z);
	se = masa_eval_source_rho_e<double>(x,y,z);

	// AD source terms
	s2u = evaluate_q(xyz,1);
	s2v = evaluate_q(xyz,2);
	s2w = evaluate_q(xyz,3);
	s2p = evaluate_q(xyz,4);
	s2e = evaluate_q(xyz,5);

	unorm = fabs(su-s2u);
	vnorm = fabs(sv-s2v);
	wnorm = fabs(sw-s2w);
	pnorm = fabs(sp-s2p);
	enorm = fabs(se-s2e);

	double urnorm = fabs(su-s2u)/std::max(su,s2u);
	double vrnorm = fabs(sv-s2v)/std::max(sv,s2v);
	double wrnorm = fabs(sw-s2w)/std::max(sw,s2w);
	double prnorm = fabs(sp-s2p)/std::max(sp,s2p);
	double ernorm = fabs(se-s2e)/std::max(se,s2e);

	unorm_max = std::max(unorm, unorm_max);
	vnorm_max = std::max(vnorm, vnorm_max);
	wnorm_max = std::max(wnorm, wnorm_max);
	pnorm_max = std::max(pnorm, pnorm_max);
	enorm_max = std::max(enorm, enorm_max);

	urnorm_max = std::max(urnorm, urnorm_max);
	vrnorm_max = std::max(vrnorm, vrnorm_max);
	wrnorm_max = std::max(wrnorm, wrnorm_max);
	prnorm_max = std::max(prnorm, prnorm_max);
	ernorm_max = std::max(ernorm, ernorm_max);

      }
    }
  }

  // std::cout << "max error in u      : " << unorm_max << std::endl;
  // std::cout << "max error in v      : " << vnorm_max << std::endl;
  // std::cout << "max error in w      : " << wnorm_max << std::endl;
  // std::cout << "max error in density: " << pnorm_max << std::endl;
  // std::cout << "max error in energy : " << enorm_max << std::endl;

  threshcheck(urnorm_max);
  threshcheck(urnorm_max);
  threshcheck(vrnorm_max);
  threshcheck(wrnorm_max);
  threshcheck(prnorm_max);
  threshcheck(ernorm_max);

  // steady as she goes...
  return 0;

}

template <std::size_t NDIM, typename ADScalar>
double evaluate_q (const NumberVector<NDIM, ADScalar>& xyz, const int ret)
{
  using std::cos;
  using std::sqrt;

  typedef typename RawType<ADScalar>::value_type Scalar;

  const Scalar PI = std::acos(Scalar(-1));

  const Scalar L = masa_get_param("L");
  const Scalar R = masa_get_param("R");
  const Scalar k = masa_get_param("k");
  const Scalar Gamma = masa_get_param("Gamma");
  const Scalar mu = masa_get_param("mu");
  const Scalar mu_bulk = masa_get_param("mu_bulk");
  const Scalar rho_0 = masa_get_param("rho_0");
  const Scalar rho_x = masa_get_param("rho_x");
  const Scalar rho_y = masa_get_param("rho_y");
  const Scalar rho_z = masa_get_param("rho_z");
  const Scalar u_0 = masa_get_param("u_0");
  const Scalar u_x = masa_get_param("u_x");
  const Scalar u_y = masa_get_param("u_y");
  const Scalar u_z = masa_get_param("u_z");
  const Scalar v_0 = masa_get_param("v_0");
  const Scalar v_x = masa_get_param("v_x");
  const Scalar v_y = masa_get_param("v_y");
  const Scalar v_z = masa_get_param("v_z");
  const Scalar w_0 = masa_get_param("w_0");
  const Scalar w_x = masa_get_param("w_x");
  const Scalar w_y = masa_get_param("w_y");
  const Scalar w_z = masa_get_param("w_z");
  const Scalar p_0 = masa_get_param("p_0");
  const Scalar p_x = masa_get_param("p_x");
  const Scalar p_y = masa_get_param("p_y");
  const Scalar p_z = masa_get_param("p_z");
  const Scalar a_rhox = masa_get_param("a_rhox");
  const Scalar a_rhoy = masa_get_param("a_rhoy");
  const Scalar a_rhoz = masa_get_param("a_rhoz");
  const Scalar a_ux = masa_get_param("a_ux");
  const Scalar a_uy = masa_get_param("a_uy");
  const Scalar a_uz = masa_get_param("a_uz");
  const Scalar a_vx = masa_get_param("a_vx");
  const Scalar a_vy = masa_get_param("a_vy");
  const Scalar a_vz = masa_get_param("a_vz");
  const Scalar a_wx = masa_get_param("a_wx");
  const Scalar a_wy = masa_get_param("a_wy");
  const Scalar a_wz = masa_get_param("a_wz");
  const Scalar a_px = masa_get_param("a_px");
  const Scalar a_py = masa_get_param("a_py");
  const Scalar a_pz = masa_get_param("a_pz");
  const Scalar Cs = masa_get_param("Cs");
  const Scalar CI = masa_get_param("CI");
  const Scalar PrT= masa_get_param("PrT");
  const Scalar deltabar = masa_get_param("deltabar");

  const ADScalar& x = xyz[0];
  const ADScalar& y = xyz[1];
  const ADScalar& z = xyz[2];

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  // Arbitrary manufactured solution
  U[0] = u_0 + u_x * cos(a_ux * PI * x / L) * u_y * cos(a_uy * PI * y / L) * u_z * cos(a_uz * PI * z / L);
  U[1] = v_0 + v_x * cos(a_vx * PI * x / L) * v_y * cos(a_vy * PI * y / L) * v_z * cos(a_vz * PI * z / L);
  U[2] = w_0 + w_x * cos(a_wx * PI * x / L) * w_y * cos(a_wy * PI * y / L) * w_z * cos(a_wz * PI * z / L);
  ADScalar RHO = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * rho_z * cos(a_rhoz * PI * z / L);
  ADScalar P = p_0 + p_x * cos(a_px * PI * x / L) * p_y * cos(a_py * PI * y / L) * p_z * cos(a_pz * PI * z / L);

  // Temperature
  ADScalar T = P / RHO / R;

  // Perfect gas energies
  ADScalar E = 1./(Gamma-1.)*P;
  ADScalar ET = E + .5 * RHO * U.dot(U);

  // The shear strain tensor
  NumberVector<NDIM, typename ADScalar::derivatives_type> GradU = gradient(U);

  // The identity tensor I
  NumberVector<NDIM, NumberVector<NDIM, Scalar> > Identity =
    NumberVector<NDIM, Scalar>::identity();

  // Constant Smagorinsky
  NumberVector<NDIM, NumberVector<NDIM, ADScalar> > S = 0.5*(GradU + transpose(GradU));
  ADScalar Smag = sqrt(2.0 * (S[0][0]*S[0][0] + S[0][1]*S[0][1] + S[0][2]*S[0][2]
			      + S[1][0]*S[1][0] + S[1][1]*S[1][1] + S[1][2]*S[1][2]
			      + S[2][0]*S[2][0] + S[2][1]*S[2][1] + S[2][2]*S[2][2]));
  ADScalar mut = (Cs*deltabar) * (Cs*deltabar) * RHO * Smag;
  ADScalar sigmakk = CI * deltabar*deltabar * RHO * Smag * Smag;

  // The shear stress tensor
  NumberVector<NDIM, NumberVector<NDIM, ADScalar> > Tau = 2.0 * (mu + mut) * (S - 1./3.*divergence(U)*Identity) + mu_bulk * divergence(U)*Identity - 2./3. * sigmakk * Identity;

  // Temperature flux
  double Cv = R / (Gamma - 1);
  NumberVector<NDIM, ADScalar> q = -(k + Gamma * Cv * mut/PrT) * T.derivatives();

  Scalar Q_rho;
  NumberVector<NDIM, Scalar> Q_rho_u;
  Scalar Q_rho_e;

  switch(ret){

    // u
  case 1:
    Q_rho_u = raw_value(divergence(RHO*U.outerproduct(U) - Tau) + P.derivatives());
    return Q_rho_u[0];
    break;

    // v
  case 2:
    Q_rho_u = raw_value(divergence(RHO*U.outerproduct(U) - Tau) + P.derivatives());
    return Q_rho_u[1];
    break;

    // w
  case 3:
    Q_rho_u = raw_value(divergence(RHO*U.outerproduct(U) - Tau) + P.derivatives());
    return Q_rho_u[2];
    break;

    // rho
  case 4:
    Q_rho = raw_value(divergence(RHO*U));
    return Q_rho;
    break;

    // energy
  case 5:
    Q_rho_e = raw_value(divergence((ET+P)*U + q - Tau.dot(U)));
    return Q_rho_e;
    break;

  default:
    std::cout << "something is wrong!\n";
    exit(1);
  }
}

#else // HAVE_METAPHYSICL

int main(void)
{
  return 77; // Autotools code for "skip test"
}

#endif // HAVE_METAPHYSICL
