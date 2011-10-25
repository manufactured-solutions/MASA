
// private, called from exact_t
double helper_func_src(double x)
{

  double helper_func;
  helper_func = x+12.4;
  return helper_func;

}

// public, but will be called from eval_q_t
double eval_exact_u(double x)
{
  double exact_u;
  exact_u = x*x;
  return exact_u;
}

double eval_q_t(double x) const
{
  double Q_T;
  Q_T = helper_func_an(x)* A_x * A_x * k_0 * std::cos(A_x * x)*eval_q_u(x);
  return Q_T;
}
