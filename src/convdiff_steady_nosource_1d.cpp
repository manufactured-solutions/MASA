// -*-c++-*-
//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include <masa_internal.h>

#include <ad_masa.h>

typedef ShadowNumber<double, long double> RawScalar;
const unsigned int NDIM = 1;
typedef DualNumber<RawScalar, NumberArray<NDIM, RawScalar> > FirstDerivType;
typedef DualNumber<FirstDerivType, NumberArray<NDIM, FirstDerivType> > SecondDerivType;
typedef SecondDerivType ADType;

using namespace MASA;

template <typename Scalar>
MASA::convdiff_steady_nosource_1d<Scalar>::convdiff_steady_nosource_1d()
{
  this->mmsname = "convdiff_steady_nosource_1d";
  this->dimension = 1;

  this->register_var("a_ux",&a_ux);
  this->register_var("a_cx",&a_cx);
  this->register_var("nu",&nu);
  this->register_var("L",&L);

  this->init_var();

} // done with constructor

template <typename Scalar>
int MASA::convdiff_steady_nosource_1d<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("a_ux",2.0);
  err += this->set_var("a_cx",2.0);
  err += this->set_var("nu",1.0);
  err += this->set_var("L",1.0);

  return err;

} // done with init_var

// ----------------------------------------
// Source Terms
// ----------------------------------------

// public static method, that can be called from eval_q_t
template <typename Scalar>
Scalar MASA::convdiff_steady_nosource_1d<Scalar>::eval_q_c(Scalar x, Scalar y) const
{
  Scalar q_c;
  q_c = PI * a_cx * std::cos(PI * a_cx * x / L) * std::cos(PI * a_ux * x / L) / L - PI * a_ux * std::sin(PI * a_cx * x / L) * std::sin(PI * a_ux * x / L) / L + std::pow(PI, 2) * std::pow(a_cx, 2) * nu * std::sin(PI * a_cx * x / L) / std::pow(L, 2);
  return q_c;
}



// ----------------------------------------
// Analytical Terms
// ----------------------------------------

// public method
template <typename Scalar>
Scalar MASA::convdiff_steady_nosource_1d<Scalar>::eval_exact_c(Scalar x)
{
  Scalar exact_c;
  exact_c = std::sin(a_cx * PI * x / L);
  return exact_c;
}

// public method
template <typename Scalar>
Scalar MASA::convdiff_steady_nosource_1d<Scalar>::eval_exact_u(Scalar x)
{
  Scalar exact_u;
  exact_u = std::cos(a_ux * PI * x / L);
  return exact_u;
}



// ----------------------------------------
// Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::convdiff_steady_nosource_1d);



//---------------------------------------------------------
// AUTOMASA
// Generated on: 2013-04-23 14:32:21
//---------------------------------------------------------
