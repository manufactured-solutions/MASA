
// example of a public method called from eval_exact_t
double eval_exact_u(double x, double y)
{
  double exact_u;
  exact_u = u_0 + u_x * std::cos(a_ux * PI * x / L) * u_y * std::cos(a_uy * PI * y / L);
  return exact_u;
}

// public method
double eval_exact_v(double x, double y)
{
  double exact_v;
  exact_v = v_0 + v_x * std::cos(a_vx * PI * x / L) * v_y * std::cos(a_vy * PI * y / L);
  return exact_v;
}

// public method
double eval_exact_p(double x, double y)
{
  ADScalar P = p_0 + p_x * std::cos(a_px * PI * x / L) * p_y * std::cos(a_py * PI * y / L);
  return P;
}

// public method
double eval_exact_rho(double x, double y)
{
  ADScalar RHO = rho_0 + rho_x * std::cos(a_rhox * PI * x / L) * rho_y * std::cos(a_rhoy * PI * y / L);
  return RHO;
}

