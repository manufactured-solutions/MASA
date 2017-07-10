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

double my_beta;
double my_gamma;
double my_delta;
double my_omega;
double kx;
double ky;
double kz;

double a;
double b;
double c;
double d;
double e;
double nu;
double grav;
double gradRho;
double diff;

// typedef double RawScalar;
typedef ShadowNumber<double, long double> RawScalar;

const unsigned int NDIM = 3;

typedef DualNumber<RawScalar, NumberVector<NDIM, RawScalar> > FirstDerivType;
typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;

typedef SecondDerivType ADType;
// typedef FirstDerivType ADType;

template <std::size_t NDIM, typename RawScalar>
double evaluate_qu (const NumberVector<NDIM, RawScalar>& xyz);

template <std::size_t NDIM, typename RawScalar>
double evaluate_qv (const NumberVector<NDIM, RawScalar>& xyz);

template <std::size_t NDIM, typename RawScalar>
double evaluate_qw (const NumberVector<NDIM, RawScalar>& xyz);

template <std::size_t NDIM, typename RawScalar>
double evaluate_qr (const NumberVector<NDIM, RawScalar>& xyz);

using namespace MASA;

int main(void)
{
  int err = 0;
  int N   = 4; // mesh pts. in x and y and z
  double su,sv,sw,sr,s2u,s2v,s2w,s2r;
  double unorm, vnorm, wnorm, rnorm;
  double unorm_max, vnorm_max, wnorm_max, rnorm_max;
  double urnorm_max = 0., vrnorm_max = 0., wrnorm_max = 0., rrnorm_max = 0.;

  unorm_max = 0;
  vnorm_max = 0;
  wnorm_max = 0;
  rnorm_max = 0;

  // initialize the problem in MASA
  err += masa_init("test_bouss","navierstokes_3d_incompbouss_homogeneous");
  err += masa_sanity_check();

  my_beta = masa_get_param("beta");
  my_gamma = masa_get_param("gamma");
  my_delta = masa_get_param("delta");
  my_omega = masa_get_param("omega");
  kx = masa_get_param("kx");
  ky = masa_get_param("ky");
  kz = masa_get_param("kz");
  a = masa_get_param("a");
  b = masa_get_param("b");
  c = masa_get_param("c");
  d = masa_get_param("d");
  e = masa_get_param("e");
  nu = masa_get_param("nu");
  grav = masa_get_param("grav");
  gradRho = masa_get_param("gradRho");
  diff = masa_get_param("diff");

  // we first set up the DualNumbers that correspond to independent
  // variables, spatial coordinates x and y and z

  NumberVector<NDIM, RawScalar> xyz;

  // the input argument xyz is another NumberVector 
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
              sr  = masa_eval_source_rho<double>(i*h,j*h,k*h);

	      // AD source terms
	      s2u = evaluate_qu(xyz);
	      s2v = evaluate_qv(xyz);
              s2w = evaluate_qw(xyz);
              s2r = evaluate_qr(xyz);

	      unorm = fabs(su-s2u);	  
	      vnorm = fabs(sv-s2v);
              wnorm = fabs(sw-s2w);
              rnorm = fabs(sr-s2r);

	      double urnorm = fabs(su-s2u)/std::max(su,s2u);	  
	      double vrnorm = fabs(sv-s2v)/std::max(sv,s2v);
              double wrnorm = fabs(sw-s2w)/std::max(sw,s2w);
              double rrnorm = fabs(sr-s2r)/std::max(sr,s2r);

	      unorm_max = std::max(unorm, unorm_max);
	      vnorm_max = std::max(vnorm, vnorm_max);
              wnorm_max = std::max(wnorm, wnorm_max);
              rnorm_max = std::max(rnorm, rnorm_max);

	      urnorm_max = std::max(urnorm, urnorm_max);
	      vrnorm_max = std::max(vrnorm, vrnorm_max);
              wrnorm_max = std::max(wrnorm, wrnorm_max);
              rrnorm_max = std::max(rrnorm, rrnorm_max);

	    }
	}
    }//done with 'k' index

  // steady as she goes...
  return 0;

}

template <typename Scalar>
Scalar helper_f(Scalar x)
{
  Scalar func;
  func = 1/(my_beta+std::sin(kx*x));
  return func;
}

template <typename Scalar>
Scalar helper_frho(Scalar x)
{
  Scalar func;
  func = 1/(my_gamma+std::sin(kx*x));
  return func;
}

template <typename Scalar>
Scalar helper_g(Scalar y)
{
  Scalar func;
  func = 1/(my_delta+std::sin(ky*y));
  return func;
}

template <typename Scalar>
Scalar helper_grho(Scalar y)
{
  Scalar func;
  func = 1/(my_omega+std::sin(ky*y));
  return func;
}

template <typename Scalar>
Scalar helper_h(Scalar z)
{
  Scalar func;
  func = 1/(my_gamma+std::sin(kz*z));
  return func;
}

template <typename Scalar>
Scalar helper_hrho(Scalar z)
{
  Scalar func;
  func = 1/(my_delta+std::sin(kz*z));
  return func;
}


// Note: ADScalar needs to be a FirstDerivType or better since we have
// first derivatives here.  Adding diffusion will require a
// SecondDerivType or better

template <std::size_t NDIM, typename RawScalar>
double evaluate_qu (const NumberVector<NDIM, RawScalar>& xyz)
{
  typedef DualNumber<RawScalar, NumberVector<NDIM, RawScalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberVector<NDIM, SecondDerivType> > ThirdDerivType;
  typedef ThirdDerivType ADScalar;

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  ADScalar x = ADScalar(xyz[0],NumberVectorUnitVector<NDIM, 0, RawScalar>::value());
  ADScalar y = ADScalar(xyz[1],NumberVectorUnitVector<NDIM, 1, RawScalar>::value());
  ADScalar z = ADScalar(xyz[2],NumberVectorUnitVector<NDIM, 2, RawScalar>::value());

  // Arbitrary manufactured solutions
  U[0]          = a * helper_f(x)                  * helper_g(y).derivatives()[1]     * helper_h(z).derivatives()[2];
  U[1]          = b * helper_f(x).derivatives()[0] * helper_g(y)                      * helper_h(z).derivatives()[2];
  U[2]          = c * helper_f(x).derivatives()[0] * helper_g(y).derivatives()[1]     * helper_h(z);
  ADScalar P    = d * helper_f(x)                  * helper_g(y)                      * helper_h(z);
  ADScalar RHO  = e * helper_frho(x)               * helper_grho(y).derivatives()[1]  * helper_hrho(z);

  // Treat gravity as a vector
  NumberVector<NDIM, ADScalar> G;

  G[0] = 0.0;
  G[1] = -grav;
  G[2] = 0.0;

  // NS equation residuals
  NumberVector<NDIM, RawScalar> Q_u = 
    raw_value(

	      // convective term
	      - divergence(U.outerproduct(U))

	      // pressure
	      - P.derivatives()

              // Boussinesq term
              + RHO*G

	      // dissipation
	      + nu * divergence(gradient(U)));

  return raw_value(Q_u[0]);
}

template <std::size_t NDIM, typename RawScalar>
double evaluate_qv (const NumberVector<NDIM, RawScalar>& xyz)
{
  typedef DualNumber<RawScalar, NumberVector<NDIM, RawScalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberVector<NDIM, SecondDerivType> > ThirdDerivType;
  typedef ThirdDerivType ADScalar;

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  ADScalar x = ADScalar(xyz[0],NumberVectorUnitVector<NDIM, 0, RawScalar>::value());
  ADScalar y = ADScalar(xyz[1],NumberVectorUnitVector<NDIM, 1, RawScalar>::value());
  ADScalar z = ADScalar(xyz[2],NumberVectorUnitVector<NDIM, 2, RawScalar>::value());

  // Arbitrary manufactured solutions
  U[0]          = a * helper_f(x)                  * helper_g(y).derivatives()[1]     * helper_h(z).derivatives()[2];
  U[1]          = b * helper_f(x).derivatives()[0] * helper_g(y)                      * helper_h(z).derivatives()[2];
  U[2]          = c * helper_f(x).derivatives()[0] * helper_g(y).derivatives()[1]     * helper_h(z);
  ADScalar P    = d * helper_f(x)                  * helper_g(y)                      * helper_h(z);
  ADScalar RHO  = e * helper_frho(x)               * helper_grho(y).derivatives()[1]  * helper_hrho(z);
  
  // Treat gravity as a vector
  NumberVector<NDIM, ADScalar> G;

  G[0] = 0.0;
  G[1] = -grav;
  G[2] = 0.0;

  // NS equation residuals
  NumberVector<NDIM, RawScalar> Q_u =
    raw_value(

              // convective term
              - divergence(U.outerproduct(U))

              // pressure
              - P.derivatives()

              // Boussinesq term
              + RHO*G

              // dissipation
              + nu * divergence(gradient(U)));

  return raw_value(Q_u[1]);
}

template <std::size_t NDIM, typename RawScalar>
double evaluate_qw (const NumberVector<NDIM, RawScalar>& xyz)
{
  typedef DualNumber<RawScalar, NumberVector<NDIM, RawScalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberVector<NDIM, SecondDerivType> > ThirdDerivType;
  typedef ThirdDerivType ADScalar;

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  ADScalar x = ADScalar(xyz[0],NumberVectorUnitVector<NDIM, 0, RawScalar>::value());
  ADScalar y = ADScalar(xyz[1],NumberVectorUnitVector<NDIM, 1, RawScalar>::value());
  ADScalar z = ADScalar(xyz[2],NumberVectorUnitVector<NDIM, 2, RawScalar>::value());

  // Arbitrary manufactured solutions
  U[0]          = a * helper_f(x)                  * helper_g(y).derivatives()[1]     * helper_h(z).derivatives()[2];
  U[1]          = b * helper_f(x).derivatives()[0] * helper_g(y)                      * helper_h(z).derivatives()[2];
  U[2]          = c * helper_f(x).derivatives()[0] * helper_g(y).derivatives()[1]     * helper_h(z);
  ADScalar P    = d * helper_f(x)                  * helper_g(y)                      * helper_h(z);
  ADScalar RHO  = e * helper_frho(x)               * helper_grho(y).derivatives()[1]  * helper_hrho(z);
  
  // Treat gravity as a vector
  NumberVector<NDIM, ADScalar> G;

  G[0] = 0.0;
  G[1] = -grav;
  G[2] = 0.0;

  // NS equation residuals
  NumberVector<NDIM, RawScalar> Q_u =
    raw_value(

              // convective term
              - divergence(U.outerproduct(U))

              // pressure
              - P.derivatives()

              // Boussinesq term
              + RHO*G

              // dissipation
              + nu * divergence(gradient(U)));

  return raw_value(Q_u[2]);
}

template <std::size_t NDIM, typename RawScalar>
double evaluate_qr (const NumberVector<NDIM, RawScalar>& xyz)
{
  typedef DualNumber<RawScalar, NumberVector<NDIM, RawScalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberVector<NDIM, SecondDerivType> > ThirdDerivType;
  typedef ThirdDerivType ADScalar;

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  ADScalar x = ADScalar(xyz[0],NumberVectorUnitVector<NDIM, 0, RawScalar>::value());
  ADScalar y = ADScalar(xyz[1],NumberVectorUnitVector<NDIM, 1, RawScalar>::value());
  ADScalar z = ADScalar(xyz[2],NumberVectorUnitVector<NDIM, 2, RawScalar>::value());

  // Arbitrary manufactured solutions
  U[0]          = a * helper_f(x)                  * helper_g(y).derivatives()[1]     * helper_h(z).derivatives()[2];
  U[1]          = b * helper_f(x).derivatives()[0] * helper_g(y)                      * helper_h(z).derivatives()[2];
  U[2]          = c * helper_f(x).derivatives()[0] * helper_g(y).derivatives()[1]     * helper_h(z);
  ADScalar P    = d * helper_f(x)                  * helper_g(y)                      * helper_h(z);
  ADScalar RHO  = e * helper_frho(x)               * helper_grho(y).derivatives()[1]  * helper_hrho(z);
  
  // NS equation residuals
  double Q_rho =
    raw_value(

              // convective term
              - divergence(RHO*U)

              // Boussinesq term
              - gradRho*U[1]

              // dissipation
              + diff * divergence(gradient(RHO)));

  return Q_rho;
}

#else // HAVE_METAPHYSICL

int main(void)
{
  return 77; // Autotools code for "skip test"
}

#endif // HAVE_METAPHYSICL
