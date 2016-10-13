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
// ad_euler.h: test the euler equation is roughly identical to 
//             the source terms generated from maple
//
// $Id: $
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------


#include <tests.h>
#include <iostream>
#include "ad_masa.h"

typedef double RawScalar;

template <std::size_t NDIM, typename Scalar>
double evaluate_q (const NumberVector<NDIM, Scalar>& xyz, const int);

using namespace MASA;

int main(void)
{
  int err = 0;
  int N   = 10; // mesh pts. in x and y
  double su,sv,s2u,s2v,sp,se,s2e,s2p;
  double pnorm, unorm, vnorm, enorm;
  double pnorm_max, unorm_max, vnorm_max, enorm_max;
  double prnorm_max = 0., urnorm_max = 0., vrnorm_max = 0., ernorm_max = 0.;

  unorm_max = 0;
  vnorm_max = 0;
  pnorm_max = 0;
  enorm_max = 0;

  const unsigned int NDIM = 2;

  const RawScalar xvecinit[] = {1., 0.};
  const RawScalar yvecinit[] = {0., 1.};

  const NumberVector<NDIM, RawScalar> xvec(xvecinit);
  const NumberVector<NDIM, RawScalar> yvec(yvecinit);

  typedef DualNumber<RawScalar, NumberVector<NDIM, RawScalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;

  typedef SecondDerivType ADType;
  // typedef FirstDerivType ADType;

  // initialize the problem in MASA
  err += masa_init("euler-maple","euler_2d");
  
  // call the sanity check routine
  // (tests that all variables have been initialized)
  err += masa_sanity_check();
  //err += masa_printid<double>();

  // we first set up the DualNumber::derivatives() terms.  
  // When main() says "xy[0] = ADType(1., xvec);", that's saying "x = 1, and 
  // the gradient of f(x,y)=x is the constant vector xvec={1,0}"  
  // Likewise "xy[1] = ADType(1., yvec);" means "y = 1, and the gradient of f(x,y)=y 
  // is the constant vector yvec={0,1}" 
  NumberVector<NDIM, ADType> xy;
  xy[0] = ADType(1., xvec);
  xy[1] = ADType(1., yvec);

  // the input argument xyz is another NumberVector 
  // a vector just like Q_rho_u, a spatial location rather 
  // than a vector-valued forcing function.
  double h = 1.0/N;
  for (int i=0; i != N+1; ++i)
    {
      xy[0] = ADType(i*h,xvec);
      for (int j=0; j != N+1; ++j)
	{
	  xy[1] = ADType(j*h,yvec);
	  // evaluate masa source terms
	  su  = masa_eval_source_rho_u<double>(i*h,j*h);
	  sv  = masa_eval_source_rho_v<double>(i*h,j*h);
	  sp  = masa_eval_source_rho  <double>(i*h,j*h);
	  se  = masa_eval_source_rho_e<double>(i*h,j*h);

	  // AD source terms
	  s2u = evaluate_q(xy,1);
	  s2v = evaluate_q(xy,2);
	  s2p = evaluate_q(xy,3);
	  s2e = evaluate_q(xy,4);

	  unorm = fabs(su-s2u);	  
	  vnorm = fabs(sv-s2v);
	  pnorm = fabs(sp-s2p);	  
	  enorm = fabs(se-s2e);

	  double urnorm = fabs(su-s2u)/std::max(su,s2u);	  
	  double vrnorm = fabs(sv-s2v)/std::max(sv,s2v);
	  double prnorm = fabs(sp-s2p)/std::max(sp,s2p);	  
	  double ernorm = fabs(se-s2e)/std::max(se,s2e);

          unorm_max = std::max(unorm, unorm_max);
          vnorm_max = std::max(vnorm, vnorm_max);
          pnorm_max = std::max(pnorm, pnorm_max);
          enorm_max = std::max(enorm, enorm_max);

          urnorm_max = std::max(urnorm, urnorm_max);
          vrnorm_max = std::max(vrnorm, vrnorm_max);
          prnorm_max = std::max(prnorm, prnorm_max);
          ernorm_max = std::max(ernorm, ernorm_max);

	}
    }
 

  threshcheck(urnorm_max);
  threshcheck(vrnorm_max);
  threshcheck(prnorm_max);
  threshcheck(ernorm_max);

  // steady as she goes...
  return 0;

}

// Note: ADScalar needs to be a FirstDerivType or better since we have
// first derivatives here.  Adding diffusion will require a
// SecondDerivType or better

template <std::size_t NDIM, typename ADScalar>
double evaluate_q (const NumberVector<NDIM, ADScalar>& xyz, const int ret)
{
  typedef typename RawType<ADScalar>::value_type Scalar;

  const Scalar PI = std::acos(Scalar(-1));

  const Scalar k = 1.38;
  const Scalar u_0 = 200.23;
  const Scalar u_x = 1.1;
  const Scalar u_y = 1.08;
  const Scalar v_0 = 1.2;
  const Scalar v_x = 1.6;
  const Scalar v_y = .47;
  const Scalar rho_0 = 100.02;
  const Scalar rho_x = 2.22;
  const Scalar rho_y = 0.8;
  const Scalar p_0 = 150.2;
  const Scalar p_x = .91;
  const Scalar p_y = .623;
  const Scalar a_px = .165;
  const Scalar a_py = .612;
  const Scalar a_rhox = 1.0;
  const Scalar a_rhoy = 1.0;
  const Scalar a_ux = .1987;
  const Scalar a_uy = 1.189;
  const Scalar a_vx = 1.91;
  const Scalar a_vy = 1.0;
  const Scalar Gamma = 1.01;
  const Scalar mu = .918;
  const Scalar L = 3.02;

  const ADScalar& x = xyz[0];
  const ADScalar& y = xyz[1];

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  // Arbitrary manufactured solution
  U[0] = u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L);
  U[1] = v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L);
  ADScalar RHO = rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L);
  ADScalar P = p_0 + p_x * std::cos(a_px * PI * x / L) + p_y * std::sin(a_py * PI * y / L);

  // Perfect gas energies
  ADScalar E = 1./(Gamma-1.)*P/RHO;
  ADScalar ET = E + .5 * U.dot(U);

  // Euler equation residuals
  Scalar Q_rho = raw_value(divergence(RHO*U));
  NumberVector<NDIM, Scalar> Q_rho_u = raw_value(divergence(RHO*U.outerproduct(U)) + P.derivatives());

  // energy equation
  Scalar Q_rho_e = raw_value(divergence((RHO*ET+P)*U));
  // Scalar Q_rho_e = raw_value(divergence((RHO*U*ET)+(P*U)));

  switch(ret)
    {

      // u
    case 1: 
      return Q_rho_u[0];
      break;

      // v
    case 2:
      return Q_rho_u[1];
      break;

      // rho
    case 3:
      return Q_rho;
      break;

      // energy
    case 4:
      return Q_rho_e;
      break;

    default:
      std::cout << "something is wrong!\n";
      exit;
    }
  return 0;
}
