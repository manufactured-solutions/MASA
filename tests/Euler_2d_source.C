#include <math.h>

const double PI = acos(-1);

double SourceQ_e (
  double x,
  double y,
  double u_0,
  double u_x,
  double u_y,
  double v_0,
  double v_x,
  double v_y,
  double rho_0,
  double rho_x,
  double rho_y,
  double p_0,
  double p_x,
  double p_y,
  double a_px,
  double a_py,
  double a_rhox,
  double a_rhoy,
  double a_ux,
  double a_uy,
  double a_vx,
  double a_vy,
  double Gamma,
  double mu,
  double L)
{
  double Q_e;
  Q_e = -Gamma * (u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0) * a_px * p_x * PI * sin(a_px * PI * x / L) / (Gamma - 0.1e1) / L + Gamma * (v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0) * a_py * p_y * PI * cos(a_py * PI * y / L) / (Gamma - 0.1e1) / L + a_rhox * PI * rho_x * cos(a_rhox * PI * x / L) * (u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0) * (pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0, 0.2e1) + pow(u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1)) / L / 0.2e1 - a_rhoy * PI * rho_y * sin(a_rhoy * PI * y / L) * (v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0) * (pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0, 0.2e1) + pow(u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1)) / L / 0.2e1 + (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_0) * a_ux * PI * u_x * cos(x * a_ux * PI / L) * Gamma / (Gamma - 0.1e1) / L + (pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0, 0.2e1) + 0.3e1 * pow(u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * a_ux * PI * u_x * cos(x * a_ux * PI / L) / L / 0.2e1 + (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_0) * a_vy * PI * v_y * cos(y * a_vy * PI / L) * Gamma / (Gamma - 0.1e1) / L + (0.3e1 * pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0, 0.2e1) + pow(u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * a_vy * PI * v_y * cos(y * a_vy * PI / L) / L / 0.2e1 - (v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * (u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0) * PI * a_uy * u_y * sin(a_uy * PI * y / L) / L - (v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * (u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0) * PI * a_vx * v_x * sin(a_vx * PI * x / L) / L;
  return(Q_e);
}

double SourceQ_u ( // 23 variables
  double x,
  double y,
  double u_0,
  double u_x,
  double u_y,
  double v_0,
  double v_x,
  double v_y,
  double rho_0,
  double rho_x,
  double rho_y,
  double p_0,
  double p_x,
  double p_y,
  double a_px,
  double a_py,
  double a_rhox,
  double a_rhoy,
  double a_ux,
  double a_uy,
  double a_vx,
  double a_vy,
  double L)
{
  double Q_u;
  Q_u = -p_x * sin(a_px * PI * x / L) * a_px * PI / L + rho_x * cos(a_rhox * PI * x / L) * pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L), 0.2e1) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_rhoy * PI / L + 0.2e1 * u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_ux * PI / L - u_y * sin(a_uy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_uy * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_vy * PI / L;
  return(Q_u);
}

double SourceQ_v (
  double x,
  double y,
  double u_0,
  double u_x,
  double u_y,
  double v_0,
  double v_x,
  double v_y,
  double rho_0,
  double rho_x,
  double rho_y,
  double p_0,
  double p_x,
  double p_y,
  double a_px,
  double a_py,
  double a_rhox,
  double a_rhoy,
  double a_ux,
  double a_uy,
  double a_vx,
  double a_vy,
  double L)
{
  double Q_v;
  Q_v = p_y * cos(a_py * PI * y / L) * a_py * PI / L + rho_x * cos(a_rhox * PI * x / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_rhox * PI / L - sin(a_rhoy * PI * y / L) * rho_y * pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L), 0.2e1) * a_rhoy * PI / L + cos(a_ux * PI * x / L) * u_x * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_ux * PI / L - sin(a_vx * PI * x / L) * v_x * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_vx * PI / L + 0.2e1 * cos(a_vy * PI * y / L) * v_y * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_vy * PI / L;
  return(Q_v);
}

double SourceQ_rho(
  double x,
  double y,
  double u_0,
  double u_x,
  double u_y,
  double v_0,
  double v_x,
  double v_y,
  double rho_0,
  double rho_x,
  double rho_y,
  double p_0,
  double p_x,
  double p_y,
  double a_px,
  double a_py,
  double a_rhox,
  double a_rhoy,
  double a_ux,
  double a_uy,
  double a_vx,
  double a_vy,
  double L)
{
  double Q_rho;
  Q_rho = rho_x * cos(a_rhox * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_rhoy * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * a_ux * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * a_vy * PI / L;
  return(Q_rho);
}
