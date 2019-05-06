// -*-c++-*-
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

#include <masa_internal.h>

#ifdef HAVE_METAPHYSICL

#include <ad_masa.h>

typedef ShadowNumber<double, long double> RawScalar;
const unsigned int NDIM = 3;
typedef DualNumber<RawScalar, NumberVector<NDIM, RawScalar> > FirstDerivType;
typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
typedef SecondDerivType ADType;

using namespace MASA;

template <typename Scalar>
MASA::ad_cns_3d_les_sph<Scalar>::ad_cns_3d_les_sph()
{
  this->mmsname = "ad_cns_3d_les_sph";
  this->dimension = 3;

  this->register_var("R",&R);
  this->register_var("k",&k);
  this->register_var("u_0",&u_0);
  this->register_var("v_0",&v_0);
  this->register_var("w_0",&w_0);
  this->register_var("u_r",&u_r);
  this->register_var("v_r",&v_r);
  this->register_var("w_r",&w_r);
  this->register_var("rho_0",&rho_0);
  this->register_var("rho_x",&rho_x);
  this->register_var("rho_y",&rho_y);
  this->register_var("rho_z",&rho_z);
  this->register_var("p_0",&p_0);
  this->register_var("p_r",&p_r);
  this->register_var("a_pr",&a_pr);
  this->register_var("a_rhox",&a_rhox);
  this->register_var("a_rhoy",&a_rhoy);
  this->register_var("a_rhoz",&a_rhoz);
  this->register_var("a_ur",&a_ur);
  this->register_var("a_vr",&a_vr);
  this->register_var("a_wr",&a_wr);
  this->register_var("Gamma",&Gamma);
  this->register_var("mu",&mu);
  this->register_var("mu_bulk",&mu_bulk);
  this->register_var("L",&L);
  this->register_var("Cs",&Cs);
  this->register_var("CI",&CI);
  this->register_var("PrT",&PrT);
  this->register_var("deltabar",&deltabar);

  this->init_var();

} // done with constructor

template <typename Scalar>
int MASA::ad_cns_3d_les_sph<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("R",1.01);
  err += this->set_var("k",1.38);
  err += this->set_var("u_0",10.23);
  err += this->set_var("v_0",10.23);
  err += this->set_var("w_0",10.23);
  err += this->set_var("u_r",10.23);
  err += this->set_var("v_r",12);
  err += this->set_var("w_r",12);
  err += this->set_var("rho_0",10.02);
  err += this->set_var("rho_x",7.2);
  err += this->set_var("rho_y",9.8);
  err += this->set_var("rho_z",9.8);
  err += this->set_var("p_0",10.2);
  err += this->set_var("p_r",.91);
  err += this->set_var("a_pr",.165);
  err += this->set_var("a_rhox",.627);
  err += this->set_var("a_rhoy",.828);
  err += this->set_var("a_rhoz",.828);
  err += this->set_var("a_ur",.1987);
  err += this->set_var("a_vr",1.91);
  err += this->set_var("a_wr",1.91);
  err += this->set_var("Gamma",1.01);
  err += this->set_var("mu",.918);
  err += this->set_var("mu_bulk",.3);
  err += this->set_var("L",3.02);
  err += this->set_var("Cs",0.16);
  err += this->set_var("CI",0.09);
  err += this->set_var("PrT",0.7);
  err += this->set_var("deltabar",1.0);

  return err;

} // done with init_var

// ----------------------------------------
// Source Terms
// ----------------------------------------

// public static method, that can be called from eval_q_t
template <typename Scalar>
Scalar MASA::ad_cns_3d_les_sph<Scalar>::eval_q_rho_u(Scalar x1, Scalar y1, Scalar z1)
{
  using std::cos;
  using std::sqrt;

  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef SecondDerivType ADScalar;

  const ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  const ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  const ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  // Arbitrary manufactured solution
  const ADScalar r = sqrt(x * x + y * y + z * z);
  U[0] = u_0 + u_r * cos(a_ur * PI * r);
  U[1] = v_0 + v_r * cos(a_vr * PI * r);
  U[2] = w_0 + w_r * cos(a_wr * PI * r);
  ADScalar RHO = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * rho_z * cos(a_rhoz * PI * z / L);
  ADScalar P = p_0 + p_r * cos(a_pr * PI * r);

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

  // Momentum equation
  NumberVector<NDIM, Scalar> Q_rho_u =
    raw_value(divergence(RHO*U.outerproduct(U) - Tau) + P.derivatives());

  return Q_rho_u[0];
}

// public, static method
template <typename Scalar>
Scalar MASA::ad_cns_3d_les_sph<Scalar>::eval_q_rho_v(Scalar x1, Scalar y1, Scalar z1)
{
  using std::cos;
  using std::sqrt;

  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef SecondDerivType ADScalar;

  const ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  const ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  const ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  // Arbitrary manufactured solution
  const ADScalar r = sqrt(x * x + y * y + z * z);
  U[0] = u_0 + u_r * cos(a_ur * PI * r);
  U[1] = v_0 + v_r * cos(a_vr * PI * r);
  U[2] = w_0 + w_r * cos(a_wr * PI * r);
  ADScalar RHO = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * rho_z * cos(a_rhoz * PI * z / L);
  ADScalar P = p_0 + p_r * cos(a_pr * PI * r);

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

  // Momentum equation
  NumberVector<NDIM, Scalar> Q_rho_u =
    raw_value(divergence(RHO*U.outerproduct(U) - Tau) + P.derivatives());

  return Q_rho_u[1];

}

// public, static method
template <typename Scalar>
Scalar MASA::ad_cns_3d_les_sph<Scalar>::eval_q_rho_w(Scalar x1, Scalar y1, Scalar z1)
{
  using std::cos;
  using std::sqrt;

  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef SecondDerivType ADScalar;

  const ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  const ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  const ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  // Arbitrary manufactured solution
  const ADScalar r = sqrt(x * x + y * y + z * z);
  U[0] = u_0 + u_r * cos(a_ur * PI * r);
  U[1] = v_0 + v_r * cos(a_vr * PI * r);
  U[2] = w_0 + w_r * cos(a_wr * PI * r);
  ADScalar RHO = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * rho_z * cos(a_rhoz * PI * z / L);
  ADScalar P = p_0 + p_r * cos(a_pr * PI * r);

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

  // Momentum equation
  NumberVector<NDIM, Scalar> Q_rho_u =
    raw_value(divergence(RHO*U.outerproduct(U) - Tau) + P.derivatives());

  return Q_rho_u[2];

}

// public, static method
template <typename Scalar>
Scalar MASA::ad_cns_3d_les_sph<Scalar>::eval_q_rho_e(Scalar x1, Scalar y1, Scalar z1)
{
  using std::cos;
  using std::sqrt;

  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef SecondDerivType ADScalar;

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  const ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  const ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  const ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  // Arbitrary manufactured solution
  const ADScalar r = sqrt(x * x + y * y + z * z);
  U[0] = u_0 + u_r * cos(a_ur * PI * r);
  U[1] = v_0 + v_r * cos(a_vr * PI * r);
  U[2] = w_0 + w_r * cos(a_wr * PI * r);
  ADScalar RHO = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * rho_z * cos(a_rhoz * PI * z / L);
  ADScalar P = p_0 + p_r * cos(a_pr * PI * r);

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

  // Energy equation
  Scalar Q_rho_e = raw_value(divergence((ET+P)*U + q - Tau.dot(U)));

  return Q_rho_e;
}


// public, static method
template <typename Scalar>
Scalar MASA::ad_cns_3d_les_sph<Scalar>::eval_q_rho(Scalar x1, Scalar y1, Scalar z1)
{
  using std::cos;
  using std::sqrt;

  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef SecondDerivType ADScalar;

  const ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  const ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  const ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  // Arbitrary manufactured solution
  const ADScalar r = sqrt(x * x + y * y + z * z);
  U[0] = u_0 + u_r * cos(a_ur * PI * r);
  U[1] = v_0 + v_r * cos(a_vr * PI * r);
  U[2] = w_0 + w_r * cos(a_wr * PI * r);
  ADScalar RHO = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * rho_z * cos(a_rhoz * PI * z / L);
  ADScalar P = p_0 + p_r * cos(a_pr * PI * r);

  // Continuity equation
  Scalar Q_rho = raw_value(divergence(RHO*U));

  return Q_rho;
}


// ----------------------------------------
// Analytical Terms
// ----------------------------------------

// example of a public method called from eval_exact_t
template <typename Scalar>
Scalar MASA::ad_cns_3d_les_sph<Scalar>::eval_exact_u(Scalar x, Scalar y, Scalar z)
{
  using std::cos;

  Scalar exact_u;
  const Scalar r = sqrt(x * x + y * y + z * z);
  exact_u = u_0 + u_r * cos(a_ur * PI * r);
  return exact_u;
}

// public method
template <typename Scalar>
Scalar MASA::ad_cns_3d_les_sph<Scalar>::eval_exact_v(Scalar x, Scalar y, Scalar z)
{
  using std::cos;

  Scalar exact_v;
  const Scalar r = sqrt(x * x + y * y + z * z);
  exact_v = v_0 + v_r * cos(a_vr * PI * r);
  return exact_v;
}

// public method
template <typename Scalar>
Scalar MASA::ad_cns_3d_les_sph<Scalar>::eval_exact_w(Scalar x, Scalar y, Scalar z)
{
  using std::cos;

  Scalar exact_w;
  const Scalar r = sqrt(x * x + y * y + z * z);
  exact_w = w_0 + w_r * cos(a_wr * PI * r);
  return exact_w;
}

// public method
template <typename Scalar>
Scalar MASA::ad_cns_3d_les_sph<Scalar>::eval_exact_p(Scalar x, Scalar y, Scalar z)
{
  using std::cos;

  const Scalar r = sqrt(x * x + y * y + z * z);
  Scalar P = p_0 + p_r * cos(a_pr * PI * r);
  return P;
}

// public method
template <typename Scalar>
Scalar MASA::ad_cns_3d_les_sph<Scalar>::eval_exact_rho(Scalar x, Scalar y, Scalar z)
{
  using std::cos;

  Scalar RHO = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * rho_z * cos(a_rhoz * PI * z / L);
  return RHO;
}



// ----------------------------------------
// Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::ad_cns_3d_les_sph);


#endif // HAVE_METAPHYSICL
