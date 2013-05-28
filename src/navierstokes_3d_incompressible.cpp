// -*-c++-*-
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

#include <masa_internal.h>

#include <ad_masa.h>



// Private methods declarations
template <typename Scalar, typename Scalar2>
Scalar helper_f(Scalar2 beta, Scalar2 kx, Scalar x);

template <typename Scalar>
Scalar helper_g(Scalar y);
  
template <typename Scalar>
Scalar helper_gt(Scalar y);

template <typename Scalar, typename Scalar2>
Scalar helper_h(Scalar2 gamma, Scalar2 kz, Scalar z);



typedef ShadowNumber<double, long double> RawScalar;
const unsigned int NDIM = 3;
typedef DualNumber<RawScalar, NumberArray<NDIM, RawScalar> > FirstDerivType;
typedef DualNumber<FirstDerivType, NumberArray<NDIM, FirstDerivType> > SecondDerivType;
typedef SecondDerivType ADType;

using namespace MASA;

template <typename Scalar>
MASA::navierstokes_3d_incompressible<Scalar>::navierstokes_3d_incompressible()
{
  this->mmsname = "navierstokes_3d_incompressible";
  this->dimension = 3;

  this->register_var("a",&a);
  this->register_var("b",&b);
  this->register_var("c",&c);
  this->register_var("d",&d);
  this->register_var("beta",&beta);
  this->register_var("gamma",&gamma);
  this->register_var("nu",&nu);
  this->register_var("kx",&kx);
  this->register_var("kz",&kz);

  this->init_var();

} // done with constructor

template <typename Scalar>
int MASA::navierstokes_3d_incompressible<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("a",1.05);
  err += this->set_var("b",2.15);
  err += this->set_var("c",-3.2);
  err += this->set_var("d",10.1);
  err += this->set_var("beta",2.2);
  err += this->set_var("gamma",2.4);
  err += this->set_var("nu",.02);
  err += this->set_var("kx",1);
  err += this->set_var("kz",1);  

  return err;

} // done with init_var

// ----------------------------------------
// Source Terms
// ----------------------------------------

// public static method, that can be called from eval_q_t
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible<Scalar>::eval_q_u(Scalar x1, Scalar y1, Scalar z1)
{
  typedef DualNumber<Scalar, NumberArray<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberArray<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberArray<NDIM, SecondDerivType> > ThirdDerivType;
  typedef ThirdDerivType ADScalar;

  // Treat velocity as a vector
  NumberArray<NDIM, ADScalar> U;

  ADScalar x = ADScalar(x1,NumberArrayUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberArrayUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberArrayUnitVector<NDIM, 2, Scalar>::value());

  // Arbitrary manufactured solutions
  U[0]       = a * helper_f(beta,kx,x)                  + helper_g(y).derivatives()[1] + helper_h(gamma,kz,z).derivatives()[2];
  U[1]       = b * helper_f(beta,kx,x).derivatives()[0] + helper_g(y)                  + helper_h(gamma,kz,z).derivatives()[2];
  U[2]       = c * helper_f(beta,kx,x).derivatives()[0] + helper_g(y).derivatives()[1] + helper_h(gamma,kz,z);
  ADScalar P = d * helper_f(beta,kx,x)                  + helper_gt(y)                 + helper_h(gamma,kz,z);

  // NS equation residuals
  NumberArray<NDIM, Scalar> Q_rho_u = 
    raw_value(

	      // convective term
	      divergence(U.outerproduct(U))

	      // pressure
	      - P.derivatives()

	      // dissipation
	      + nu * divergence(gradient(U)));

  return Q_rho_u[0];
}


// public static method, that can be called from eval_q_t
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible<Scalar>::eval_q_v(Scalar x1, Scalar y1, Scalar z1)
{
  typedef DualNumber<Scalar, NumberArray<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberArray<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberArray<NDIM, SecondDerivType> > ThirdDerivType;
  typedef ThirdDerivType ADScalar;

  // Treat velocity as a vector
  NumberArray<NDIM, ADScalar> U;

  ADScalar x = ADScalar(x1,NumberArrayUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberArrayUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberArrayUnitVector<NDIM, 2, Scalar>::value());

  // Arbitrary manufactured solutions
  U[0]       = a * helper_f(beta,kx,x)                  + helper_g(y).derivatives()[1] + helper_h(gamma,kz,z).derivatives()[2];
  U[1]       = b * helper_f(beta,kx,x).derivatives()[0] + helper_g(y)                  + helper_h(gamma,kz,z).derivatives()[2];
  U[2]       = c * helper_f(beta,kx,x).derivatives()[0] + helper_g(y).derivatives()[1] + helper_h(gamma,kz,z);
  ADScalar P = d * helper_f(beta,kx,x)                  + helper_gt(y)                 + helper_h(gamma,kz,z);

  // NS equation residuals
  NumberArray<NDIM, Scalar> Q_rho_u = 
    raw_value(

	      // convective term
	      divergence(U.outerproduct(U))

	      // pressure
	      - P.derivatives()

	      // dissipation
	      + nu * divergence(gradient(U)));

  return Q_rho_u[1];
}


// public static method, that can be called from eval_q_t
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible<Scalar>::eval_q_w(Scalar x1, Scalar y1, Scalar z1)
{
  typedef DualNumber<Scalar, NumberArray<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberArray<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberArray<NDIM, SecondDerivType> > ThirdDerivType;
  typedef ThirdDerivType ADScalar;

  // Treat velocity as a vector
  NumberArray<NDIM, ADScalar> U;

  ADScalar x = ADScalar(x1,NumberArrayUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberArrayUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberArrayUnitVector<NDIM, 2, Scalar>::value());

  // Arbitrary manufactured solutions
  U[0]       = a * helper_f(beta,kx,x)                  + helper_g(y).derivatives()[1] + helper_h(gamma,kz,z).derivatives()[2];
  U[1]       = b * helper_f(beta,kx,x).derivatives()[0] + helper_g(y)                  + helper_h(gamma,kz,z).derivatives()[2];
  U[2]       = c * helper_f(beta,kx,x).derivatives()[0] + helper_g(y).derivatives()[1] + helper_h(gamma,kz,z);
  ADScalar P = d * helper_f(beta,kx,x)                  + helper_gt(y)                 + helper_h(gamma,kz,z);

  // NS equation residuals
  NumberArray<NDIM, Scalar> Q_rho_u = 
    raw_value(

	      // convective term
	      divergence(U.outerproduct(U))

	      // pressure
	      - P.derivatives()

	      // dissipation
	      + nu * divergence(gradient(U)));

  return Q_rho_u[2];
}

// ----------------------------------------
// Analytical Terms
// ----------------------------------------

// example of a private method, called from exact_t
template <typename Scalar, typename Scalar2>
Scalar helper_f(Scalar2 beta, Scalar2 kx, Scalar x)
//Scalar MASA::navierstokes_3d_incompressible<Scalar>::helper_f(Scalar x)
{
  Scalar func;
  func = 1/(beta+std::sin(kx*x));
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

template <typename Scalar, typename Scalar2>
//Scalar MASA::navierstokes_3d_incompressible<Scalar>::helper_h(Scalar z)
Scalar helper_h(Scalar2 gamma, Scalar2 kz, Scalar z)
{
  Scalar func;
  func = 1/(gamma+std::sin(kz*z));
  return func;
}

//
// main functions
// 

// example of a public method called from eval_exact_t
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible<Scalar>::eval_exact_u(Scalar x, Scalar y1, Scalar z1)
{
  typedef DualNumber<Scalar, Scalar> OneDDerivType;
  OneDDerivType y = OneDDerivType(y1,1);
  OneDDerivType z = OneDDerivType(z1,1);

  Scalar exact_u;
  exact_u =   a *  helper_f(beta,kx,x) + helper_g(y).derivatives() +  helper_h(gamma,kz,z).derivatives();
  return exact_u;
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible<Scalar>::eval_exact_v(Scalar x1, Scalar y, Scalar z1)
{
  typedef DualNumber<Scalar, Scalar> OneDDerivType;
  OneDDerivType x = OneDDerivType(x1,1);
  OneDDerivType z = OneDDerivType(z1,1);

  Scalar exact_v;
  exact_v = b * helper_f(beta,kx,x).derivatives() +  helper_g(y) + helper_h(gamma,kz,z).derivatives();
  return exact_v;
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible<Scalar>::eval_exact_w(Scalar x1, Scalar y1, Scalar z)
{
  typedef DualNumber<Scalar, Scalar> OneDDerivType;
  OneDDerivType x = OneDDerivType(x1,1);
  OneDDerivType y = OneDDerivType(y1,1);

  Scalar exact_w;
  exact_w = c * helper_f(beta,kx,x).derivatives() + helper_g(y).derivatives() +  helper_h(gamma,kz,z);
  return exact_w;
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible<Scalar>::eval_exact_p(Scalar x, Scalar y, Scalar z)
{
  Scalar P = d *  helper_f(beta,kx,x) + helper_gt(y) +  helper_h(gamma,kz,z);
  return P;
}


// ----------------------------------------
// Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::navierstokes_3d_incompressible);



//---------------------------------------------------------
// AUTOMASA
// Generated on: 2013-05-06 12:39:39
//---------------------------------------------------------
