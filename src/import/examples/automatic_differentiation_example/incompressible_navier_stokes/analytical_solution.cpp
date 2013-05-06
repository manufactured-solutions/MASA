
// example of a private method, called from exact_t
double helper_f(double x)
{
  double func;
  func = 1/(beta+std::sin(kx*x);
  return func;
}

double helper_g(double y)
{
  double func;
  func = (1-y*y)*(1-y*y)*gt(y);
  return func;
}
  
double helper_gt(double y)
{

}

double helper_h(double z)
{
  double func;
  func = 1/(gamma+std::sin(kz*z);
  return func;
}

//
// derivatives of above functions
//
double helper_fp(double x)
{

}

double helper_gp(double y)
{

}

double helper_hp(double z)
{

}

//
// main functions
// 

// example of a public method called from eval_exact_t
double eval_exact_u(double x, double y, double z)
{
  double exact_u;
  exact_u =   a *  helper_f(x) + helper_gp(y) +  helper_h(z);
  return exact_u;
}

// public method
double eval_exact_v(double x, double y, double z)
{
  double exact_v;
  exact_v = b * helper_fp(x) +  helper_g(y) + helper_hp(z);
  return exact_v;
}

// public method
double eval_exact_w(double x, double y, double z)
{
  double exact_w;
  exact_w = c * helper_fp(x) + helper_gp(y) +  helper_h(z);
  return exact_w;
}

// public method
double eval_exact_p(double x, double y, double z)
{
  double P = d *  helper_f(x) + helper_gt(y) +  helper_h(z);
  return P;
}
