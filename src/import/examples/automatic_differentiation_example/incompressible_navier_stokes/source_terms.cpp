
// u component of velocity source term
double eval_q_u(double x, double y, double z)
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
  ADScalar P = d * helper_f(x)                  + helper_gt(y)                 + helper_h(z);

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
double eval_q_v(double x, double y, double z)
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
  ADScalar P = d * helper_f(x)                  + helper_gt(y)                 + helper_h(z);

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
double eval_q_w(double x, double y, double z)
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
  ADScalar P = d * helper_f(x)                  + helper_gt(y)                 + helper_h(z);

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

