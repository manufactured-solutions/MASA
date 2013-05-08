// -*-c++-*-
//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include <masa_internal.h>

#include <ad_masa.h>

typedef ShadowNumber<double, long double> RawScalar;
const unsigned int NDIM = 3;
typedef DualNumber<RawScalar, NumberArray<NDIM, RawScalar> > FirstDerivType;
typedef DualNumber<FirstDerivType, NumberArray<NDIM, FirstDerivType> > SecondDerivType;
typedef SecondDerivType ADType;

using namespace MASA;

double ad_beta_hom;
double ad_gamma_hom;
double ad_delta_hom;
double ad_kx_hom;
double ad_kz_hom;
double ad_ky_hom;

template <typename Scalar>
MASA::navierstokes_3d_incompressible_homogeneous<Scalar>::navierstokes_3d_incompressible_homogeneous()
{
  this->mmsname = "navierstokes_3d_incompressible_homogeneous";
  this->dimension = 3;

  this->register_var("a",&a);
  this->register_var("b",&b);
  this->register_var("c",&c);
  this->register_var("d",&d);
  this->register_var("beta",&beta);
  this->register_var("gamma",&gamma);
  this->register_var("delta",&delta);
  this->register_var("nu",&nu);
  this->register_var("kx",&kx);
  this->register_var("kz",&kz);
  this->register_var("ky",&ky);

  this->init_var();

} // done with constructor

template <typename Scalar>
int MASA::navierstokes_3d_incompressible_homogeneous<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("a",1.05);
  err += this->set_var("b",2.15);
  err += this->set_var("c",-3.2);
  err += this->set_var("d",10.1);
  err += this->set_var("beta",2.2);
  err += this->set_var("gamma",2.4);
  err += this->set_var("delta",2.0);
  err += this->set_var("nu",.02);
  err += this->set_var("kx",1);
  err += this->set_var("kz",1);
  err += this->set_var("ky",1);

  return err;

} // done with init_var

// ----------------------------------------
// Source Terms
// ----------------------------------------

// u component of velocity source term
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible_homogeneous<Scalar>::eval_q_u(Scalar x1, Scalar y1, Scalar z1)
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
  U[0]       = a * helper_f(x)                  + helper_g(y).derivatives()[1] + helper_h(z).derivatives()[2];
  U[1]       = b * helper_f(x).derivatives()[0] + helper_g(y)                  + helper_h(z).derivatives()[2];
  U[2]       = c * helper_f(x).derivatives()[0] + helper_g(y).derivatives()[1] + helper_h(z);
  ADScalar P = d * helper_f(x)                  + helper_g(y)                  + helper_h(z);

  // NS equation residuals
  NumberArray<NDIM, Scalar> Q_u = 
    raw_value(

	      // convective term
	      divergence(U.outerproduct(U))

	      // pressure
	      - P.derivatives()

	      // dissipation
	      + nu * divergence(gradient(U)));

  return Q_u[1];

}

// v component of velocity source term
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible_homogeneous<Scalar>::eval_q_v(Scalar x1, Scalar y1, Scalar z1)
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
  U[0]       = a * helper_f(x)                  + helper_g(y).derivatives()[1] + helper_h(z).derivatives()[2];
  U[1]       = b * helper_f(x).derivatives()[0] + helper_g(y)                  + helper_h(z).derivatives()[2];
  U[2]       = c * helper_f(x).derivatives()[0] + helper_g(y).derivatives()[1] + helper_h(z);
  ADScalar P = d * helper_f(x)                  + helper_g(y)                  + helper_h(z);

  // NS equation residuals
  NumberArray<NDIM, Scalar> Q_u = 
    raw_value(

	      // convective term
	      divergence(U.outerproduct(U))

	      // pressure
	      - P.derivatives()

	      // dissipation
	      + nu * divergence(gradient(U)));

  return Q_u[1];

}

// w component of velocity source term
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible_homogeneous<Scalar>::eval_q_w(Scalar x1, Scalar y1, Scalar z1)
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
  U[0]       = a * helper_f(x)                  + helper_g(y).derivatives()[1] + helper_h(z).derivatives()[2];
  U[1]       = b * helper_f(x).derivatives()[0] + helper_g(y)                  + helper_h(z).derivatives()[2];
  U[2]       = c * helper_f(x).derivatives()[0] + helper_g(y).derivatives()[1] + helper_h(z);
  ADScalar P = d * helper_f(x)                  + helper_g(y)                  + helper_h(z);

  // NS equation residuals
  NumberArray<NDIM, Scalar> Q_u = 
    raw_value(

	      // convective term
	      divergence(U.outerproduct(U))

	      // pressure
	      - P.derivatives()

	      // dissipation
	      + nu * divergence(gradient(U)));

  return Q_u[2];

}



// ----------------------------------------
// Analytical Terms
// ----------------------------------------

// helper functions
template <typename Scalar>
Scalar helper_f(Scalar x)
{
  Scalar func;
  func = 1/(ad_beta_hom+std::sin(ad_kx_hom*x));
  return func;
}

template <typename Scalar>
Scalar helper_g(Scalar y)
{
  Scalar func;
  func = 1/(ad_delta_hom+std::sin(ad_ky_hom*y));
  return func;
}
  
template <typename Scalar>
Scalar helper_h(Scalar z)
{
  Scalar func;
  func = 1/(ad_gamma_hom+std::sin(ad_kz_hom*z));
  return func;
}

//
// main functions
// 

// example of a public method called from eval_exact_t
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible_homogeneous<Scalar>::eval_exact_u(Scalar x, Scalar y1, Scalar z)
{
  typedef DualNumber<Scalar, Scalar> OneDDerivType;
  OneDDerivType y = OneDDerivType(y1,1);
 
  Scalar exact_u;
  exact_u =   a *  helper_f(x) + helper_g(y).derivatives() +  helper_h(z);
  return exact_u;
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible_homogeneous<Scalar>::eval_exact_v(Scalar x1, Scalar y, Scalar z1)
{
  typedef DualNumber<Scalar, Scalar> OneDDerivType;
  OneDDerivType x = OneDDerivType(x1,1);
  OneDDerivType z = OneDDerivType(z1,1);

  Scalar exact_v;
  exact_v = b * helper_f(x).derivatives() +  helper_g(y) + helper_h(z).derivatives();
  return exact_v;
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible_homogeneous<Scalar>::eval_exact_w(Scalar x1, Scalar y1, Scalar z)
{
  typedef DualNumber<Scalar, Scalar> OneDDerivType;
  OneDDerivType x = OneDDerivType(x1,1);
  OneDDerivType y = OneDDerivType(y1,1);

  Scalar exact_w;
  exact_w = c * helper_f(x).derivatives() + helper_g(y).derivatives() +  helper_h(z);
  return exact_w;
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible_homogeneous<Scalar>::eval_exact_p(Scalar x, Scalar y, Scalar z)
{

  Scalar P = d *  helper_f(x) + helper_g(y) +  helper_h(z);
  return P;
}


// ----------------------------------------
// Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::navierstokes_3d_incompressible_homogeneous);



//---------------------------------------------------------
// AUTOMASA
// Generated on: 2013-05-08 11:32:28
//---------------------------------------------------------
