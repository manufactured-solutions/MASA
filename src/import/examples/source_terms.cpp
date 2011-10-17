
double MASA::heateq_1d_steady_const<double>::eval_q_t(double x)
{
  double Q_T;
  Q_T = A_x * A_x * k_0 * std::cos(A_x * x);
  return Q_T;
}
