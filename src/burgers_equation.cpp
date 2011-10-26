// -*-c++-*-
//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include <masa_internal.h>
using namespace MASA;

template <typename Scalar>
MASA::burgers_equation<Scalar>::burgers_equation()
{
  this->mmsname = "burgers_equation";
  this->dimension = 2;

  this->register_var("nu",&nu);
  this->register_var("u_0",&u_0);
  this->register_var("u_x",&u_x);
  this->register_var("u_y",&u_y);
  this->register_var("u_t",&u_t);
  this->register_var("v_0",&v_0);
  this->register_var("v_x",&v_x);
  this->register_var("v_y",&v_y);
  this->register_var("v_t",&v_t);
  this->register_var("a_ux",&a_ux);
  this->register_var("a_uy",&a_uy);
  this->register_var("a_ut",&a_ut);
  this->register_var("a_vx",&a_vx);
  this->register_var("a_vy",&a_vy);
  this->register_var("a_vt",&a_vt);
  this->register_var("L",&L);

  this->init_var();

} // done with constructor

template <typename Scalar>
int MASA::burgers_equation<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("nu",12);
  err += this->set_var("u_0",12);
  err += this->set_var("u_x",12);
  err += this->set_var("u_y",12);
  err += this->set_var("u_t",12);
  err += this->set_var("v_0",12);
  err += this->set_var("v_x",12);
  err += this->set_var("v_y",12);
  err += this->set_var("v_t",12);
  err += this->set_var("a_ux",12);
  err += this->set_var("a_uy",12);
  err += this->set_var("a_ut",12);
  err += this->set_var("a_vx",12);
  err += this->set_var("a_vy",12);
  err += this->set_var("a_vt",12);
  err += this->set_var("L",12);

  return err;

} // done with init_var

// ----------------------------------------
// Source Terms
// ----------------------------------------
#include <math.h>

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_q_v_transient_viscous (Scalar x, Scalar y, Scalar t)
{
  Scalar Qv_tv;
  Scalar U;
  Scalar V;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Qv_tv = a_vt * PI * v_t * cos(a_vt * PI * t / L) / L - a_vx * PI * v_x * U * sin(a_vx * PI * x / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * V / L + a_vx * a_vx * PI * PI * v_x * nu * cos(a_vx * PI * x / L) * pow(L, -0.2e1) + a_vy * a_vy * PI * PI * v_y * nu * sin(a_vy * PI * y / L) * pow(L, -0.2e1);
  return(Qv_tv);
}
#include <math.h>

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_q_v_steady_viscous (Scalar x, Scalar y)
{
  Scalar Qv_sv;
  Scalar U;
  Scalar V;
  U = u_0; //+ u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0; //+ v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Qv_sv = -a_vx * PI * v_x * U * sin(a_vx * PI * x / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * V / L + a_vx * a_vx * PI * PI * v_x * nu * cos(a_vx * PI * x / L) * pow(L, -0.2e1) + a_vy * a_vy * PI * PI * v_y * nu * sin(a_vy * PI * y / L) * pow(L, -0.2e1);
  return(Qv_sv);
}
#include <math.h>

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_q_v_transient_inviscid (Scalar x, Scalar y, Scalar t)
{
  Scalar Qv_t;
  Scalar U;
  Scalar V;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Qv_t = a_vt * PI * v_t * cos(a_vt * PI * t / L) / L - a_vx * PI * v_x * U * sin(a_vx * PI * x / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * V / L;
  return(Qv_t);
}
#include <math.h>

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_q_v_steady_inviscid (Scalar x, Scalar y)
{
  Scalar Qv_s;
  Scalar U;
  Scalar V;
  U = u_0;// + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0;// + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Qv_s = -a_vx * PI * v_x * U * sin(a_vx * PI * x / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * V / L;
  return(Qv_s);
}

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_q_u_transient_viscous (Scalar x,Scalar y,Scalar t)
{
  Scalar Qu_tv;
  Scalar U;
  Scalar V;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Qu_tv = -a_ut * PI * u_t * sin(a_ut * PI * t / L) / L - a_uy * PI * u_y * V * sin(a_uy * PI * y / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * U / L + a_ux * a_ux * PI * PI * u_x * nu * sin(a_ux * PI * x / L) * pow(L, -0.2e1) + a_uy * a_uy * PI * PI * u_y * nu * cos(a_uy * PI * y / L) * pow(L, -0.2e1);
  return(Qu_tv);
}
#include <math.h>

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_q_u_steady_viscous (Scalar x, Scalar y)
{
  Scalar Qu_sv;
  Scalar U;
  Scalar V;
  U = u_0;// + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0;// + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Qu_sv = -a_uy * PI * u_y * V * sin(a_uy * PI * y / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * U / L + a_ux * a_ux * PI * PI * u_x * nu * sin(a_ux * PI * x / L) * pow(L, -0.2e1) + a_uy * a_uy * PI * PI * u_y * nu * cos(a_uy * PI * y / L) * pow(L, -0.2e1);
  return(Qu_sv);
}
#include <math.h>

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_q_u_transient_inviscid (Scalar x, Scalar y, Scalar t)
{
  Scalar Qu_t;
  Scalar U;
  Scalar V;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Qu_t = -a_ut * PI * u_t * sin(a_ut * PI * t / L) / L - a_uy * PI * u_y * V * sin(a_uy * PI * y / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * U / L;
  return(Qu_t);
}
#include <math.h>

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_q_u_steady_inviscid (Scalar x, Scalar y)
{
  Scalar Qu_s;
  Scalar U;
  Scalar V;
  U = u_0;// + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0;// + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Qu_s = -a_uy * PI * u_y * V * sin(a_uy * PI * y / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * U / L;
  return(Qu_s);
}


// ----------------------------------------
// Analytical Terms
// ----------------------------------------
// public, but will be called from eval_exact_t
template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_exact_u(Scalar x,Scalar y)
{
  Scalar exact_u;
  exact_u = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  return exact_u;
}

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_exact_v(Scalar x,Scalar y)
{
  Scalar exact_v;
  exact_v = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);
  return exact_v;
}

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_exact_u_t(Scalar x,Scalar y, Scalar t)
{
  Scalar exact_u_t;
  exact_u_t = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_t * cos(a_ut * pi * t / L);
  return exact_u_t;
}

template <typename Scalar>
Scalar MASA::burgers_equation<Scalar>::eval_exact_v_t(Scalar x,Scalar y, Scalar t)
{
  Scalar exact_v_t;
  exact_v_t = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_t * sin(a_vt * pi * t / L);
  return exact_v_t;
}


// ----------------------------------------
// Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::burgers_equation);



//---------------------------------------------------------
// AUTOMASA
// Generated on: 2011-10-26 17:35:23
//---------------------------------------------------------
