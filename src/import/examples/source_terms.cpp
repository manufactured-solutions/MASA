
double eval_q_t(double x) const
{
  double Q_T;
  Q_T = A_x * A_x * k_0 * std::cos(A_x * x);
  return Q_T;
}
