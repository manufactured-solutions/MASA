
/*
// private, called from exact_t
helper_func(double x)
{



}
*/


// public, but will be called from exact_t
double eval_exact_u(double x)
{
  double exact_u;
  exact_u = x*x;
  return exact_u;
}

double eval_exact_t(double x)
{
  double exact_t;
  exact_t = std::cos(A_x * x)*eval_exact_u(x);
  //exact_t = helper_func(x)*std::cos(A_x * x);
  return exact_t;
}

