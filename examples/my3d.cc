#include <math.h>

#define PI 3.14
#define k 3.14
#define R 3.14

double SourceQ_e (
  double x,
  double y,
  double z,
  double u_0,
  double u_x,
  double u_y,
  double u_z,
  double v_0,
  double v_x,
  double v_y,
  double v_z,
  double w_0,
  double w_x,
  double w_y,
  double w_z,
  double rho_0,
  double rho_x,
  double rho_y,
  double rho_z,
  double p_0,
  double p_x,
  double p_y,
  double p_z,
  double a_px,
  double a_py,
  double a_pz,
  double a_rhox,
  double a_rhoy,
  double a_rhoz,
  double a_ux,
  double a_uy,
  double a_uz,
  double a_vx,
  double a_vy,
  double a_vz,
  double a_wx,
  double a_wy,
  double a_wz,
  double mu,
  double Gamma,
  double L)
{
   const double unique0 = pow(L,(-2));
   const double unique1 = pow(PI,2);
   const double unique4 = pow(L,(-1));
   const double unique102 = (PI * a_rhoz * z * unique4);
   const double unique101 = sin(unique102);
   const double unique105 = (PI * a_px * x * unique4);
   const double unique104 = sin(unique105);
   const double unique107 = (PI * a_py * y * unique4);
   const double unique106 = cos(unique107);
   const double unique111 = (PI * a_pz * z * unique4);
   const double unique110 = sin(unique111);
   const double unique113 = (PI * a_uy * y * unique4);
   const double unique112 = sin(unique113);
   const double unique117 = (PI * a_wy * y * unique4);
   const double unique116 = cos(unique117);
   const double unique3 = (PI * a_ux * x * unique4);
   const double unique122 = sin(unique3);
   const double unique127 = (PI * a_uz * z * unique4);
   const double unique126 = cos(unique127);
   const double unique133 = (PI * a_wx * x * unique4);
   const double unique132 = sin(unique133);
   const double unique82 = (PI * a_vy * y * unique4);
   const double unique152 = sin(unique82);
   const double unique151 = (v_y * unique152);
   const double unique86 = (PI * a_vx * x * unique4);
   const double unique85 = cos(unique86);
   const double unique84 = (v_x * unique85);
   const double unique91 = (PI * a_vz * z * unique4);
   const double unique90 = sin(unique91);
   const double unique89 = (v_z * unique90);
   const double unique150 = (v_0 + unique151 + unique84 + unique89);
   const double unique96 = (PI * a_rhox * x * unique4);
   const double unique153 = cos(unique96);
   const double unique99 = (PI * a_rhoy * y * unique4);
   const double unique154 = sin(unique99);
   const double unique155 = cos(unique102);
   const double unique156 = sin(unique86);
   const double unique157 = cos(unique91);
   const double unique121 = (u_x * unique122);
   const double unique125 = (u_z * unique126);
   const double unique162 = cos(unique113);
   const double unique161 = (u_y * unique162);
   const double unique160 = (u_0 + unique121 + unique125 + unique161);
   const double unique159 = pow(unique160,2);
   const double unique131 = (w_x * unique132);
   const double unique167 = sin(unique117);
   const double unique166 = (w_y * unique167);
   const double unique80 = (PI * a_wz * z * unique4);
   const double unique169 = cos(unique80);
   const double unique168 = (w_z * unique169);
   const double unique165 = (w_0 + unique131 + unique166 + unique168);
   const double unique164 = pow(unique165,2);
   const double unique170 = sin(unique127);
   const double unique174 = cos(unique105);
   const double unique173 = (p_x * unique174);
   const double unique176 = sin(unique107);
   const double unique175 = (p_y * unique176);
   const double unique178 = cos(unique111);
   const double unique177 = (p_z * unique178);
   const double unique172 = (p_0 + unique173 + unique175 + unique177);
   const double unique66 = (-1 + Gamma);
   const double unique65 = pow(unique66,(-1));
   const double unique171 = (Gamma * unique172 * unique4 * unique65);
   const double unique100 = (rho_z * unique101);
   const double unique95 = sin(unique96);
   const double unique94 = (rho_x * unique95);
   const double unique98 = cos(unique99);
   const double unique97 = (rho_y * unique98);
   const double unique93 = (rho_0 + unique100 + unique94 + unique97);
   const double unique179 = pow(unique93,(-3));
   const double unique180 = cos(unique133);
   const double unique181 = pow(unique93,(-1));
   const double unique183 = pow(unique150,2);
   const double unique182 = (unique159 + unique164 + unique183);
   const double unique2 = cos(unique3);
   const double unique29 = pow(R,(-1));
   const double unique79 = sin(unique80);
   const double unique81 = cos(unique82);
   const double unique92 = pow(unique93,(-2));
   const double Q_e = (-(1.333333333333333e+00 * a_ux * a_wz * mu * u_x * w_z * unique0 * unique1 * unique2 * unique79) - (1.333333333333333e+00 * a_vy * a_wz * mu * v_y * w_z * unique0 * unique1 * unique79 * unique81) - (1.333333333333333e+00 * mu * v_y * (-(unique150 * unique152) + (v_y * pow(unique81,2))) * pow(a_vy,2) * unique0 * unique1) - (2 * a_px * a_rhox * k * p_x * rho_x * unique0 * unique1 * unique104 * unique153 * unique29 * unique92) - (2 * a_py * a_rhoy * k * p_y * rho_y * unique0 * unique1 * unique106 * unique154 * unique29 * unique92) - (2 * a_pz * a_rhoz * k * p_z * rho_z * unique0 * unique1 * unique110 * unique155 * unique29 * unique92) - (2 * a_uy * a_vx * mu * u_y * v_x * unique0 * unique1 * unique112 * unique156) - (2 * a_vz * a_wy * mu * v_z * w_y * unique0 * unique1 * unique116 * unique157) - (5.000000000000000e-01 * PI * a_rhoy * rho_y * unique150 * unique154 * unique182 * unique4) - (Gamma * PI * a_px * p_x * unique104 * unique160 * unique4 * unique65) - (Gamma * PI * a_pz * p_z * unique110 * unique165 * unique4 * unique65) - (PI * a_uy * u_y * unique112 * unique150 * unique160 * unique4 * unique93) - (PI * a_uz * u_z * unique160 * unique165 * unique170 * unique4 * unique93) - (a_vx * a_vy * v_x * v_y * ((5.000000000000000e-01 * ((3 * unique183) + unique159 + unique164) * unique4 * unique93) + unique171) * unique1 * unique150 * unique156 * unique160 * unique4 * unique81 * unique93) - (k * rho_x * ((2 * rho_x * pow(unique153,2)) + (unique93 * unique95)) * pow(a_rhox,2) * unique0 * unique1 * unique172 * unique179 * unique29) - (k * rho_y * ((2 * rho_y * pow(unique154,2)) + (unique93 * unique98)) * pow(a_rhoy,2) * unique0 * unique1 * unique172 * unique179 * unique29) - (k * rho_z * ((2 * rho_z * pow(unique155,2)) + (unique101 * unique93)) * pow(a_rhoz,2) * unique0 * unique1 * unique172 * unique179 * unique29) - (mu * v_x * (-(unique150 * unique85) + (v_x * pow(unique156,2))) * pow(a_vx,2) * unique0 * unique1) - (mu * v_z * (-(unique150 * unique90) + (v_z * pow(unique157,2))) * pow(a_vz,2) * unique0 * unique1) + (1.333333333333333e+00 * a_ux * a_vy * mu * u_x * v_y * unique0 * unique1 * unique2 * unique81) + (1.333333333333333e+00 * mu * u_x * (-(u_x * pow(unique2,2)) + (unique122 * unique160)) * pow(a_ux,2) * unique0 * unique1) + (1.333333333333333e+00 * mu * w_z * (-(w_z * pow(unique79,2)) + (unique165 * unique169)) * pow(a_wz,2) * unique0 * unique1) + (2 * a_uz * a_wx * mu * u_z * w_x * unique0 * unique1 * unique170 * unique180) + (5.000000000000000e-01 * PI * a_rhox * rho_x * unique153 * unique160 * unique182 * unique4) + (5.000000000000000e-01 * a_rhoz * a_ux * rho_z * u_x * ((5.000000000000000e-01 * ((3 * unique159) + unique164 + unique183) * unique4 * unique93) + unique171) * unique1 * unique155 * unique165 * unique182 * unique2 * unique4) + (Gamma * PI * a_py * p_y * unique106 * unique150 * unique4 * unique65) + (PI * a_vz * v_z * unique150 * unique157 * unique165 * unique4 * unique93) + (PI * a_wx * w_x * unique160 * unique165 * unique180 * unique4 * unique93) + (a_wy * a_wz * w_y * w_z * (-(5.000000000000000e-01 * ((3 * unique164) + unique159 + unique183) * unique4 * unique93) - unique171) * unique1 * unique116 * unique150 * unique165 * unique4 * unique79 * unique93) + (k * p_x * pow(a_px,2) * unique0 * unique1 * unique174 * unique181 * unique29) + (k * p_y * pow(a_py,2) * unique0 * unique1 * unique176 * unique181 * unique29) + (k * p_z * pow(a_pz,2) * unique0 * unique1 * unique178 * unique181 * unique29) + (mu * u_y * (-(u_y * pow(unique112,2)) + (unique160 * unique162)) * pow(a_uy,2) * unique0 * unique1) + (mu * u_z * (-(u_z * pow(unique170,2)) + (unique126 * unique160)) * pow(a_uz,2) * unique0 * unique1) + (mu * w_x * (-(w_x * pow(unique180,2)) + (unique132 * unique165)) * pow(a_wx,2) * unique0 * unique1) + (mu * w_y * (-(w_y * pow(unique116,2)) + (unique165 * unique167)) * pow(a_wy,2) * unique0 * unique1));
  return(Q_e);
}
