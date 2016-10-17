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
#include <tests.h>
#include "ad_masa.h"

double my_beta;
double my_gamma;
double kx;
double kz;

double a;
double b;
double c;
double d;
double nu;

// typedef double RawScalar;
typedef ShadowNumber<double, long double> RawScalar;

const unsigned int NDIM = 3;

typedef DualNumber<RawScalar, NumberArray<NDIM, RawScalar> > FirstDerivType;
typedef DualNumber<FirstDerivType, NumberArray<NDIM, FirstDerivType> > SecondDerivType;

typedef SecondDerivType ADType;
// typedef FirstDerivType ADType;

template <std::size_t NDIM, typename RawScalar>
double evaluate_q (const NumberArray<NDIM, RawScalar>& xyz);

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

  // initialize the problem in MASA
  err += masa_init("ad_incom","navierstokes_3d_incompressible");
  err += masa_sanity_check();

  my_beta = masa_get_param("beta");
  my_gamma = masa_get_param("gamma");
  kx = masa_get_param("kx");
  kz = masa_get_param("kz");
  a = masa_get_param("a");
  b = masa_get_param("b");
  c = masa_get_param("c");
  d = masa_get_param("d");
  nu = masa_get_param("nu");

  // we first set up the DualNumbers that correspond to independent
  // variables, spatial coordinates x and y and z

  NumberArray<NDIM, RawScalar> xyz;

  // the input argument xyz is another NumberArray 
  // a vector just like Q_rho_u, a spatial location rather 
  // than a vector-valued forcing function.
  double h = 1.0/N;
  for (int k=0; k != N+1; ++k)
    {
      xyz[2] = k*h;

      for (int i=0; i != N+1; ++i)
	{
	  //
	  xyz[0] = i*h;
	  
	  // Under the hood:
	  // xy[0] = ADType(FirstDerivType(i*h, xvec), xvec);

	  for (int j=0; j != N+1; ++j)
	    {
	      xyz[1] = j*h;

	      // evaluate masa source terms
	      su  = masa_eval_source_u<double>(i*h,j*h,k*h);
	      sv  = masa_eval_source_v<double>(i*h,j*h,k*h);
	      sw  = masa_eval_source_w<double>(i*h,j*h,k*h);

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
  func = 1/(my_beta+std::sin(kx*x));
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
  func = 1/(my_gamma+std::sin(kz*z));
  return func;
}

// Note: ADScalar needs to be a FirstDerivType or better since we have
// first derivatives here.  Adding diffusion will require a
// SecondDerivType or better

template <std::size_t NDIM, typename RawScalar>
double evaluate_q (const NumberArray<NDIM, RawScalar>& xyz)
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
	      - divergence(U.outerproduct(U))

	      // pressure
	      - P.derivatives()

	      // dissipation
	      + nu * divergence(gradient(U)));

  return raw_value(Q_rho_u[0]);
}
