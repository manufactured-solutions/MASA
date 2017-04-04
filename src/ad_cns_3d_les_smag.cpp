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
MASA::ad_cns_3d_les_smag<Scalar>::ad_cns_3d_les_smag()
{
  this->mmsname = "ad_cns_3d_les_smag";
  this->dimension = 3;

  this->register_var("R",&R);
  this->register_var("k",&k);
  this->register_var("u_0",&u_0);
  this->register_var("u_x",&u_x);
  this->register_var("u_y",&u_y);
  this->register_var("u_z",&u_z);
  this->register_var("v_0",&v_0);
  this->register_var("v_x",&v_x);
  this->register_var("v_y",&v_y);
  this->register_var("v_z",&v_z);
  this->register_var("w_0",&w_0);
  this->register_var("w_x",&w_x);
  this->register_var("w_y",&w_y);
  this->register_var("w_z",&w_z);
  this->register_var("rho_0",&rho_0);
  this->register_var("rho_x",&rho_x);
  this->register_var("rho_y",&rho_y);
  this->register_var("rho_z",&rho_z);
  this->register_var("p_0",&p_0);
  this->register_var("p_x",&p_x);
  this->register_var("p_y",&p_y);
  this->register_var("p_z",&p_z);
  this->register_var("a_px",&a_px);
  this->register_var("a_py",&a_py);
  this->register_var("a_pz",&a_pz);
  this->register_var("a_rhox",&a_rhox);
  this->register_var("a_rhoy",&a_rhoy);
  this->register_var("a_rhoz",&a_rhoz);
  this->register_var("a_ux",&a_ux);
  this->register_var("a_uy",&a_uy);
  this->register_var("a_uz",&a_uz);
  this->register_var("a_vx",&a_vx);
  this->register_var("a_vy",&a_vy);
  this->register_var("a_vz",&a_vz);
  this->register_var("a_wx",&a_wx);
  this->register_var("a_wy",&a_wy);
  this->register_var("a_wz",&a_wz);
  this->register_var("Gamma",&Gamma);
  this->register_var("mu",&mu);
  this->register_var("L",&L);
  this->register_var("Cs",&Cs);
  this->register_var("CI",&CI);
  this->register_var("PrT",&PrT);

  this->init_var();

} // done with constructor

template <typename Scalar>
int MASA::ad_cns_3d_les_smag<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("R",1.01);
  err += this->set_var("k",1.38);
  err += this->set_var("u_0",10.23);
  err += this->set_var("u_x",1.1);
  err += this->set_var("u_y",.08);
  err += this->set_var("u_z",.08);
  err += this->set_var("v_0",12);
  err += this->set_var("v_x",1.6);
  err += this->set_var("v_y",.67);
  err += this->set_var("v_z",.67);
  err += this->set_var("w_0",12);
  err += this->set_var("w_x",1.6);
  err += this->set_var("w_y",.67);
  err += this->set_var("w_z",.67);
  err += this->set_var("rho_0",10.02);
  err += this->set_var("rho_x",7.2);
  err += this->set_var("rho_y",9.8);
  err += this->set_var("rho_z",9.8);
  err += this->set_var("p_0",10.2);
  err += this->set_var("p_x",.91);
  err += this->set_var("p_y",.623);
  err += this->set_var("p_z",.623);
  err += this->set_var("a_px",.165);
  err += this->set_var("a_py",.612);
  err += this->set_var("a_pz",.612);
  err += this->set_var("a_rhox",.627);
  err += this->set_var("a_rhoy",.828);
  err += this->set_var("a_rhoz",.828);
  err += this->set_var("a_ux",.1987);
  err += this->set_var("a_uy",1.189);
  err += this->set_var("a_uz",1.189);
  err += this->set_var("a_vx",1.91);
  err += this->set_var("a_vy",2.901);
  err += this->set_var("a_vz",2.901);
  err += this->set_var("a_wx",1.91);
  err += this->set_var("a_wy",2.901);
  err += this->set_var("a_wz",2.901);
  err += this->set_var("Gamma",1.01);
  err += this->set_var("mu",.918);
  err += this->set_var("L",3.02);
  err += this->set_var("Cs",0.16);
  err += this->set_var("CI",0.09);
  err += this->set_var("PrT",0.7);

  return err;

} // done with init_var

// ----------------------------------------
// Source Terms
// ----------------------------------------

// public static method, that can be called from eval_q_t
template <typename Scalar>
Scalar MASA::ad_cns_3d_les_smag<Scalar>::eval_q_u(Scalar x1, Scalar y1, Scalar z1) const
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
  U[0] = u_0 + u_x * cos(a_ux * PI * x / L) * u_y * cos(a_uy * PI * y / L) * cos(a_uy * PI * z / L);
  U[1] = v_0 + v_x * cos(a_vx * PI * x / L) * v_y * cos(a_vy * PI * y / L) * cos(a_vy * PI * z / L);
  U[2] = w_0 + w_x * cos(a_wx * PI * x / L) * w_y * cos(a_wy * PI * y / L) * cos(a_wy * PI * z / L);
  ADScalar RHO = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * cos(a_rhoz * PI * z / L);
  ADScalar P = p_0 + p_x * cos(a_px * PI * x / L) * p_y * cos(a_py * PI * y / L) * cos(a_pz * PI * z / L);

  // Temperature
  ADScalar T = P / RHO / R;

  // Perfect gas energies
  ADScalar E = 1./(Gamma-1.)*P/RHO;
  ADScalar ET = E + .5 * U.dot(U);

  // The shear strain tensor
  NumberVector<NDIM, typename ADScalar::derivatives_type> GradU = gradient(U);

  // The identity tensor I
  NumberVector<NDIM, NumberVector<NDIM, Scalar> > Identity = 
    NumberVector<NDIM, Scalar>::identity();

  // Constant Smagorinsky
  NumberVector<NDIM, NumberVector<NDIM, ADScalar> > S = GradU + transpose(GradU);
  ADScalar Smag = sqrt(2.0 * (S[0,0]*S[0,0] + S[0,1]*S[0,1] + S[0,2]*S[0,2] \
			      + S[1,0]*S[1,0] + S[1,1]*S[1,1] + S[1,2]*S[1,2] \
			      + S[2,0]*S[2,0] + S[2,1]*S[2,1] + S[2,2]*S[2,2] ));
  ADScalar mut = - 2.0 * (Cs*deltabar) * (Cs*deltabar) * RHO * Smag;
  ADScalar sigmakk = 2.0 * CI * deltabar*deltabar * RHO * Smag * Smag;
    
  // The shear stress tensor
  NumberVector<NDIM, NumberVector<NDIM, ADScalar> > Tau = (mu + mut) * (GradU + transpose(GradU) - 2./3.*divergence(U)*Identity) + 1./3. * sigmakk * Identity;
    
  // Temperature flux
  NumberVector<NDIM, ADScalar> q = -(k + mut/PrT) * T.derivatives();

  // Euler equation residuals
  // Scalar Q_rho = raw_value(divergence(RHO*U));
  NumberVector<NDIM, Scalar> Q_rho_u = 
    raw_value(divergence(RHO*U.outerproduct(U) - Tau) + P.derivatives());

  return Q_rho_u[0];
}

// public, static method
template <typename Scalar>
Scalar MASA::ad_cns_3d_les_smag<Scalar>::eval_q_v(Scalar x1, Scalar y1, Scalar z1) const
{
  using std::cos;

  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef SecondDerivType ADScalar;

  const ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  const ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  const ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  // Arbitrary manufactured solution
  U[0] = u_0 + u_x * cos(a_ux * PI * x / L) * u_y * cos(a_uy * PI * y / L) * cos(a_uy * PI * z / L);
  U[1] = v_0 + v_x * cos(a_vx * PI * x / L) * v_y * cos(a_vy * PI * y / L) * cos(a_vy * PI * z / L);
  U[2] = w_0 + w_x * cos(a_wx * PI * x / L) * w_y * cos(a_wy * PI * y / L) * cos(a_wy * PI * z / L);
  ADScalar RHO = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * cos(a_rhoz * PI * z / L);
  ADScalar P = p_0 + p_x * cos(a_px * PI * x / L) * p_y * cos(a_py * PI * y / L) * cos(a_pz * PI * z / L);

  // Temperature
  ADScalar T = P / RHO / R;

 // Perfect gas energies
  ADScalar E = 1./(Gamma-1.)*P/RHO;
  ADScalar ET = E + .5 * U.dot(U);

  // The shear strain tensor
  NumberVector<NDIM, typename ADScalar::derivatives_type> GradU = gradient(U);

  // The identity tensor I
  NumberVector<NDIM, NumberVector<NDIM, Scalar> > Identity = 
    NumberVector<NDIM, Scalar>::identity();

  // Constant Smagorinsky
  NumberVector<NDIM, NumberVector<NDIM, ADScalar> > S = GradU + transpose(GradU);
  ADScalar Smag = sqrt(2.0 * (S[0,0]*S[0,0] + S[0,1]*S[0,1] + S[0,2]*S[0,2] \
			      + S[1,0]*S[1,0] + S[1,1]*S[1,1] + S[1,2]*S[1,2] \
			      + S[2,0]*S[2,0] + S[2,1]*S[2,1] + S[2,2]*S[2,2] ));
  ADScalar mut = - 2.0 * (Cs*deltabar) * (Cs*deltabar) * RHO * Smag;
  ADScalar sigmakk = 2.0 * CI * deltabar*deltabar * RHO * Smag * Smag;
    
  // The shear stress tensor
  NumberVector<NDIM, NumberVector<NDIM, ADScalar> > Tau = (mu + mut) * (GradU + transpose(GradU) - 2./3.*divergence(U)*Identity) + 1./3. * sigmakk * Identity;
    
  // Temperature flux
  NumberVector<NDIM, ADScalar> q = -(k + mut/PrT) * T.derivatives();

  // Euler equation residuals
  // Scalar Q_rho = raw_value(divergence(RHO*U));
  NumberVector<NDIM, Scalar> Q_rho_u = 
    raw_value(divergence(RHO*U.outerproduct(U) - Tau) + P.derivatives());

  return Q_rho_u[1];

}

// public, static method
template <typename Scalar>
Scalar MASA::ad_cns_3d_les_smag<Scalar>::eval_q_w(Scalar x1, Scalar y1, Scalar z1) const
{
  using std::cos;

  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef SecondDerivType ADScalar;

  const ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  const ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  const ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  // Arbitrary manufactured solution
  U[0] = u_0 + u_x * cos(a_ux * PI * x / L) * u_y * cos(a_uy * PI * y / L) * cos(a_uy * PI * z / L);
  U[1] = v_0 + v_x * cos(a_vx * PI * x / L) * v_y * cos(a_vy * PI * y / L) * cos(a_vy * PI * z / L);
  U[2] = w_0 + w_x * cos(a_wx * PI * x / L) * w_y * cos(a_wy * PI * y / L) * cos(a_wy * PI * z / L);
  ADScalar RHO = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * cos(a_rhoz * PI * z / L);
  ADScalar P = p_0 + p_x * cos(a_px * PI * x / L) * p_y * cos(a_py * PI * y / L) * cos(a_pz * PI * z / L);

  // Temperature
  ADScalar T = P / RHO / R;

 // Perfect gas energies
  ADScalar E = 1./(Gamma-1.)*P/RHO;
  ADScalar ET = E + .5 * U.dot(U);

  // The shear strain tensor
  NumberVector<NDIM, typename ADScalar::derivatives_type> GradU = gradient(U);

  // The identity tensor I
  NumberVector<NDIM, NumberVector<NDIM, Scalar> > Identity = 
    NumberVector<NDIM, Scalar>::identity();

  // Constant Smagorinsky
  NumberVector<NDIM, NumberVector<NDIM, ADScalar> > S = GradU + transpose(GradU);
  ADScalar Smag = sqrt(2.0 * (S[0,0]*S[0,0] + S[0,1]*S[0,1] + S[0,2]*S[0,2] \
			      + S[1,0]*S[1,0] + S[1,1]*S[1,1] + S[1,2]*S[1,2] \
			      + S[2,0]*S[2,0] + S[2,1]*S[2,1] + S[2,2]*S[2,2] ));
  ADScalar mut = - 2.0 * (Cs*deltabar) * (Cs*deltabar) * RHO * Smag;
  ADScalar sigmakk = 2.0 * CI * deltabar*deltabar * RHO * Smag * Smag;
    
  // The shear stress tensor
  NumberVector<NDIM, NumberVector<NDIM, ADScalar> > Tau = (mu + mut) * (GradU + transpose(GradU) - 2./3.*divergence(U)*Identity) + 1./3. * sigmakk * Identity;
    
  // Temperature flux
  NumberVector<NDIM, ADScalar> q = -(k + mut/PrT) * T.derivatives();

  // Euler equation residuals
  // Scalar Q_rho = raw_value(divergence(RHO*U));
  NumberVector<NDIM, Scalar> Q_rho_u = 
    raw_value(divergence(RHO*U.outerproduct(U) - Tau) + P.derivatives());

  return Q_rho_u[2];

}

// public, static method
template <typename Scalar>
Scalar MASA::ad_cns_3d_les_smag<Scalar>::eval_q_e(Scalar x1, Scalar y1, Scalar z1) const
{
  using std::cos;

  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef SecondDerivType ADScalar;

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  const ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  const ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  const ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  // Arbitrary manufactured solution
  U[0] = u_0 + u_x * cos(a_ux * PI * x / L) * u_y * cos(a_uy * PI * y / L) * cos(a_uy * PI * z / L);
  U[1] = v_0 + v_x * cos(a_vx * PI * x / L) * v_y * cos(a_vy * PI * y / L) * cos(a_vy * PI * z / L);
  U[2] = w_0 + w_x * cos(a_wx * PI * x / L) * w_y * cos(a_wy * PI * y / L) * cos(a_wy * PI * z / L);
  ADScalar RHO = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * cos(a_rhoz * PI * z / L);
  ADScalar P = p_0 + p_x * cos(a_px * PI * x / L) * p_y * cos(a_py * PI * y / L) * cos(a_pz * PI * z / L);

  // Temperature
  ADScalar T = P / RHO / R;

  // Perfect gas energies
  ADScalar E = 1./(Gamma-1.)*P/RHO;
  ADScalar ET = E + .5 * U.dot(U);

  // The shear strain tensor
  NumberVector<NDIM, typename ADScalar::derivatives_type> GradU = gradient(U);

  // The identity tensor I
  NumberVector<NDIM, NumberVector<NDIM, Scalar> > Identity = 
    NumberVector<NDIM, Scalar>::identity();

  // Constant Smagorinsky
  NumberVector<NDIM, NumberVector<NDIM, ADScalar> > S = GradU + transpose(GradU);
  ADScalar Smag = sqrt(2.0 * (S[0,0]*S[0,0] + S[0,1]*S[0,1] + S[0,2]*S[0,2] \
			      + S[1,0]*S[1,0] + S[1,1]*S[1,1] + S[1,2]*S[1,2] \
			      + S[2,0]*S[2,0] + S[2,1]*S[2,1] + S[2,2]*S[2,2] ));
  ADScalar mut = - 2.0 * (Cs*deltabar) * (Cs*deltabar) * RHO * Smag;
  ADScalar sigmakk = 2.0 * CI * deltabar*deltabar * RHO * Smag * Smag;
    
  // The shear stress tensor
  NumberVector<NDIM, NumberVector<NDIM, ADScalar> > Tau = (mu + mut) * (GradU + transpose(GradU) - 2./3.*divergence(U)*Identity) + 1./3. * sigmakk * Identity;
    
  // Temperature flux
  NumberVector<NDIM, ADScalar> q = -(k + mut/PrT) * T.derivatives();

  // Euler equation residuals
  // Scalar Q_rho = raw_value(divergence(RHO*U));
  // NumberVector<NDIM, Scalar> Q_rho_u = 
  //   raw_value(divergence(RHO*U.outerproduct(U) - Tau) + P.derivatives());

  // energy equation
  Scalar Q_rho_e = raw_value(divergence((RHO*ET+P)*U + q - Tau.dot(U)));

  return Q_rho_e;
}


// public, static method
template <typename Scalar>
Scalar MASA::ad_cns_3d_les_smag<Scalar>::eval_q_rho(Scalar x1, Scalar y1, Scalar z1) const
{
  using std::cos;

  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef SecondDerivType ADScalar;

  const ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  const ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  const ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  // Arbitrary manufactured solution
  U[0] = u_0 + u_x * cos(a_ux * PI * x / L) * u_y * cos(a_uy * PI * y / L) * cos(a_uy * PI * z / L);
  U[1] = v_0 + v_x * cos(a_vx * PI * x / L) * v_y * cos(a_vy * PI * y / L) * cos(a_vy * PI * z / L);
  U[2] = w_0 + w_x * cos(a_wx * PI * x / L) * w_y * cos(a_wy * PI * y / L) * cos(a_wy * PI * z / L);
  ADScalar RHO = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * cos(a_rhoz * PI * z / L);
  ADScalar P = p_0 + p_x * cos(a_px * PI * x / L) * p_y * cos(a_py * PI * y / L) * cos(a_pz * PI * z / L);

  Scalar Q_rho = raw_value(divergence(RHO*U));

  return Q_rho;
}


// ----------------------------------------
// Analytical Terms
// ----------------------------------------

// example of a public method called from eval_exact_t
template <typename Scalar>
Scalar MASA::ad_cns_3d_les_smag<Scalar>::eval_exact_u(Scalar x, Scalar y, Scalar z)
{
  using std::cos;

  Scalar exact_u;
  exact_u = u_0 + u_x * cos(a_ux * PI * x / L) * u_y * cos(a_uy * PI * y / L) * cos(a_uz * PI * z / L);
  return exact_u;
}

// public method
template <typename Scalar>
Scalar MASA::ad_cns_3d_les_smag<Scalar>::eval_exact_v(Scalar x, Scalar y, Scalar z)
{
  using std::cos;

  Scalar exact_v;
  exact_v = v_0 + v_x * cos(a_vx * PI * x / L) * v_y * cos(a_vy * PI * y / L) * cos(a_vz * PI * z / L);
  return exact_v;
}

// public method
template <typename Scalar>
Scalar MASA::ad_cns_3d_les_smag<Scalar>::eval_exact_w(Scalar x, Scalar y, Scalar z)
{
  using std::cos;

  Scalar exact_v;
  exact_v = v_0 + v_x * cos(a_vx * PI * x / L) * v_y * cos(a_vy * PI * y / L) * cos(a_wz * PI * z / L);
  return exact_v;
}

// public method
template <typename Scalar>
Scalar MASA::ad_cns_3d_les_smag<Scalar>::eval_exact_p(Scalar x, Scalar y, Scalar z)
{
  using std::cos;

  Scalar P = p_0 + p_x * cos(a_px * PI * x / L) * p_y * cos(a_py * PI * y / L) * cos(a_pz * PI * z / L);
  return P;
}

// public method
template <typename Scalar>
Scalar MASA::ad_cns_3d_les_smag<Scalar>::eval_exact_rho(Scalar x, Scalar y, Scalar z)
{
  using std::cos;

  Scalar RHO = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * cos(a_rhoz * PI * z / L);
  return RHO;
}



// ----------------------------------------
// Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::ad_cns_3d_les_smag);


#endif // HAVE_METAPHYSICL
