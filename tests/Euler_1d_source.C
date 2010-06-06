#include <math.h>

#define PI acos(-1)

double SourceQ_rho (
  double x,
  double u_0,
  double u_x,
  double rho_0,
  double rho_x,
  double p_0,
  double p_x,
  double a_px,
  double a_rhox,
  double a_ux,
  double L)
{
  double Q_rho;
  Q_rho = rho_x * cos(a_rhox * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L)) * a_rhox * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L)) * a_ux * PI / L;
  return(Q_rho);
}
#include <math.h>

double SourceQ_u (
  double x,
  double u_0,
  double u_x,
  double rho_0,
  double rho_x,
  double p_0,
  double p_x,
  double a_px,
  double a_rhox,
  double a_ux,
  double L)
{
  double Q_u;
  Q_u = -sin(a_px * PI * x / L) * a_px * PI * p_x / L + rho_x * cos(a_rhox * PI * x / L) * pow(u_0 + u_x * sin(a_ux * PI * x / L), 0.2e1) * a_rhox * PI / L + 0.2e1 * u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L)) * (u_0 + u_x * sin(a_ux * PI * x / L)) * a_ux * PI / L;
  return(Q_u);
}
#include <math.h>

double SourceQ_e (
  double x,
  double u_0,
  double u_x,
  double rho_0,
  double rho_x,
  double p_0,
  double p_x,
  double a_px,
  double a_rhox,
  double a_ux,
  double Gamma,
  double mu,
  double L)
{
  double Q_e;
  Q_e = cos(a_rhox * PI * x / L) * rho_x * pow(u_0 + u_x * sin(a_ux * PI * x / L), 0.3e1) * a_rhox * PI / L / 0.2e1 + cos(a_ux * PI * x / L) * (p_0 + p_x * cos(a_px * PI * x / L)) * a_ux * PI * u_x * Gamma / L / (Gamma - 0.1e1) - Gamma * p_x * sin(a_px * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L)) * a_px * PI / L / (Gamma - 0.1e1) + 0.3e1 / 0.2e1 * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L)) * pow(u_0 + u_x * sin(a_ux * PI * x / L), 0.2e1) * a_ux * PI * u_x / L;
  return(Q_e);
}
