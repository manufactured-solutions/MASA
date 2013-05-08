// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012 The PECOS Development Team
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
#include <tests.h>
#include "ad_masa.h"

double ad_beta;
double ad_gamma;
double ad_kx;
double ad_kz;

double a = 1.05;
double b = 2.15;
double c = -3.2;
double d = 10.1;
double nu = .02;

// typedef double RawScalar;
typedef ShadowNumber<double, long double> RawScalar;

const unsigned int NDIM = 3;

typedef DualNumber<RawScalar, NumberArray<NDIM, RawScalar> > FirstDerivType;
typedef DualNumber<FirstDerivType, NumberArray<NDIM, FirstDerivType> > SecondDerivType;

typedef SecondDerivType ADType;
// typedef FirstDerivType ADType;

template <std::size_t NDIM, typename Scalar>
double evaluate_q (const NumberArray<NDIM, Scalar>& xyz);

using namespace MASA;

int main(void)
{
  int err = 0;
  int N   = 4; // mesh pts. in x and y and z
  double su,sv,sw,s2u,s2v,sp,se,s2e,s2p;
  double pnorm, unorm, vnorm, enorm;
  double pnorm_max, unorm_max, vnorm_max, enorm_max;
  double prnorm_max = 0., urnorm_max = 0., vrnorm_max = 0., ernorm_max = 0.;

  unorm_max = 0;
  vnorm_max = 0;

  const RawScalar xvecinit[] = {1., 0., 0.};
  const RawScalar yvecinit[] = {0., 1., 0.};
  const RawScalar zvecinit[] = {0., 0., 1.};

  const NumberArray<NDIM, RawScalar> xvec(xvecinit);
  const NumberArray<NDIM, RawScalar> yvec(yvecinit);
  const NumberArray<NDIM, RawScalar> zvec(zvecinit);

  // initialize the problem in MASA
  err += masa_init("ad_incom","navierstokes_3d_incompressible");
  err += masa_sanity_check();

  // we first set up the DualNumbers that correspond to independent
  // variables, spatial coordinates x and y and z

  NumberArray<NDIM, ADType> xyz;

  // When main() says "xy[0] = ADType(1., xvec);", that's saying "x = 1, and 
  // the gradient of f(x,y)=x is the constant vector xvec={1,0}"  
  // Likewise "xy[1] = ADType(1., yvec);" means "y = 1, and the gradient of f(x,y)=y 
  xyz[0] = ADType(1., xvec);
  xyz[1] = ADType(1., yvec);
  xyz[2] = ADType(1., zvec);

  // For getting second derivatives, the way to set up a
  // twice-differentiable independent variable used to be more
  // complicated: first set up a once-differentiable variable, then
  // make sure *its* derivatives are also properly set.

  // However, if the new DualNumber constructors are working properly
  // then this is unnecessary.

  // xy[0] = ADType(FirstDerivType(1., xvec), xvec);
  // xy[1] = ADType(FirstDerivType(1., yvec), yvec);

  // the input argument xyz is another NumberArray 
  // a vector just like Q_rho_u, a spatial location rather 
  // than a vector-valued forcing function.
  double h = 1.0/N;
  for (int k=0; k != N+1; ++k)
    {
      for (int i=0; i != N+1; ++i)
	{
	  //
	  xyz[0] = ADType(i*h, xvec);
	  
	  // Under the hood:
	  // xy[0] = ADType(FirstDerivType(i*h, xvec), xvec);

	  for (int j=0; j != N+1; ++j)
	    {
	      xyz[1] = ADType(j*h, yvec);

	      // evaluate masa source terms
	      su  = masa_eval_source_u<double>(k*h,i*h,j*h);
	      sv  = masa_eval_source_v<double>(k*h,i*h,j*h);
	      sw  = masa_eval_source_w<double>(k*h,i*h,j*h);

	      // AD source terms
	      s2u = evaluate_q(xyz);
	      s2v = evaluate_q(xyz);

	      unorm = fabs(su-s2u);	  
	      vnorm = fabs(sv-s2v);

	      double urnorm = fabs(su-s2u)/std::max(su,s2u);	  
	      double vrnorm = fabs(sv-s2v)/std::max(sv,s2v);

	      unorm_max = std::max(unorm, unorm_max);
	      vnorm_max = std::max(vnorm, vnorm_max);

	      urnorm_max = std::max(urnorm, urnorm_max);
	      vrnorm_max = std::max(vrnorm, vrnorm_max);

	    }
	}
    }//done with 'k' index

  // steady as she goes...
  return 0;

}

// example of a private method, called from exact_t
template <typename Scalar>
Scalar helper_f(Scalar x)
//Scalar MASA::navierstokes_3d_incompressible<Scalar>::helper_f(Scalar x)
{
  Scalar func;
  func = 1/(ad_beta+std::sin(ad_kx*x));
  return func;
}

template <typename Scalar>
Scalar helper_g(Scalar y)
//Scalar MASA::navierstokes_3d_incompressible<Scalar>::helper_g(Scalar y)
{
  Scalar func;
  func = (1-y*y)*(1-y*y)*helper_gt(y);
  return func;
}
  
template <typename Scalar>
Scalar helper_gt(Scalar y)
//Scalar MASA::navierstokes_3d_incompressible<Scalar>::helper_gt(Scalar y)
{
  Scalar func;
  func = std::pow(y,Scalar(5.0))+std::pow(y,Scalar(4.0))+std::pow(y,Scalar(3.0))+std::pow(y,Scalar(2.0))+y;
  return func;
}

template <typename Scalar>
//Scalar MASA::navierstokes_3d_incompressible<Scalar>::helper_h(Scalar z)
Scalar helper_h(Scalar z)
{
  Scalar func;
  func = 1/(ad_gamma+std::sin(ad_kz*z));
  return func;
}

// Note: ADScalar needs to be a FirstDerivType or better since we have
// first derivatives here.  Adding diffusion will require a
// SecondDerivType or better

template <std::size_t NDIM, typename ADScalar>
double evaluate_q (const NumberArray<NDIM, ADScalar>& xyz)
{
  typedef DualNumber<RawScalar, NumberArray<NDIM, RawScalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberArray<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberArray<NDIM, SecondDerivType> > ThirdDerivType;
  typedef ThirdDerivType ADScalar;

  // Treat velocity as a vector
  NumberArray<NDIM, ADScalar> U;

  ADScalar x = ADScalar(xyz[0],NumberArrayUnitVector<NDIM, 0, RawScalar>::value());
  ADScalar y = ADScalar(xyz[1],NumberArrayUnitVector<NDIM, 1, RawScalar>::value());
  ADScalar z = ADScalar(xyz[2],NumberArrayUnitVector<NDIM, 2, RawScalar>::value());

  // Arbitrary manufactured solutions
  U[0]       = a * helper_f(x)                  + helper_g(y).derivatives()[1] + helper_h(z).derivatives()[2];
  U[1]       = b * helper_f(x).derivatives()[0] + helper_g(y)                  + helper_h(z).derivatives()[2];
  U[2]       = c * helper_f(x).derivatives()[0] + helper_g(y).derivatives()[1] + helper_h(z);
  ADScalar P = d * helper_f(x)                  + helper_gt(y)                 + helper_h(z);

  // NS equation residuals
  NumberArray<NDIM, RawScalar> Q_rho_u = 
    raw_value(

	      // convective term
	      divergence(U.outerproduct(U))

	      // pressure
	      + P.derivatives()

	      // dissipation
	      + nu * divergence(gradient(U)));

  return raw_value(Q_rho_u[0]);
}
