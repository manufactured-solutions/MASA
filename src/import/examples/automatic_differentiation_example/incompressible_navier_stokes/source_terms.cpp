
// public static method, that can be called from eval_q_t
double eval_q_u(double x, double y, double z) const
{
  typedef DualNumber<Scalar, NumberArray<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberArray<NDIM, FirstDerivType> > SecondDerivType;
  typedef SecondDerivType ADScalar;

  // Treat velocity as a vector
  NumberArray<NDIM, ADScalar> U;

  // Arbitrary manufactured solutions
  U[0]       = a *  helper_f(x) + helper_gp(y) +  helper_h(z);
  U[1]       = b * helper_fp(x) +  helper_g(y) + helper_hp(z);
  U[2]       = c * helper_fp(x) + helper_gp(y) +  helper_h(z);
  ADScalar P = d *  helper_f(x) + helper_gt(y) +  helper_h(z);

  // NS equation residuals
  NumberArray<NDIM, Scalar> Q_rho_u = 
    raw_value(

	      // convective term
	      divergence(U.outerproduct(U))

	      // pressure
	      + P.derivatives()

	      // dissipation
	      + nu * divergence(gradient(U)));

  return Q_rho_u[0];
}

// public, static method
double eval_q_v(double x, double y, double z) const
{
  typedef DualNumber<Scalar, NumberArray<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberArray<NDIM, FirstDerivType> > SecondDerivType;
  typedef SecondDerivType ADScalar;

  // Treat velocity as a vector
  NumberArray<NDIM, ADScalar> U;

  // Arbitrary manufactured solutions
  U[0]       = a *  helper_f(x) + helper_gp(y) +  helper_h(z);
  U[1]       = b * helper_fp(x) +  helper_g(y) + helper_hp(z);
  U[2]       = c * helper_fp(x) + helper_gp(y) +  helper_h(z);
  ADScalar P = d *  helper_f(x) + helper_gt(y) +  helper_h(z);

  // NS equation residuals
  NumberArray<NDIM, Scalar> Q_rho_u = 
    raw_value(

	      // convective term
	      divergence(U.outerproduct(U))

	      // pressure
	      + P.derivatives()

	      // dissipation
	      + nu * divergence(gradient(U)));

  return Q_rho_u[1];

}

// public, static method
double eval_q_w(double x, double y, double z) const
{
  typedef DualNumber<Scalar, NumberArray<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberArray<NDIM, FirstDerivType> > SecondDerivType;
  typedef SecondDerivType ADScalar;

  // Treat velocity as a vector
  NumberArray<NDIM, ADScalar> U;

  // Arbitrary manufactured solutions
  U[0]       = a *  helper_f(x) + helper_gp(y) +  helper_h(z);
  U[1]       = b * helper_fp(x) +  helper_g(y) + helper_hp(z);
  U[2]       = c * helper_fp(x) + helper_gp(y) +  helper_h(z);
  ADScalar P = d *  helper_f(x) + helper_gt(y) +  helper_h(z);

  // NS equation residuals
  NumberArray<NDIM, Scalar> Q_rho_u = 
    raw_value(

	      // convective term
	      divergence(U.outerproduct(U))

	      // pressure
	      + P.derivatives()

	      // dissipation
	      + nu * divergence(gradient(U)));

  return Q_rho_u[2];

}

