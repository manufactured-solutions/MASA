// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// $Author: 
// $Id:
//
// euler_chem.cpp: 
//          These are the MASA class member functions and constructors
//          For the Euler Equations + Chemistry
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h> 

using namespace MASA;
/* ------------------------------------------------
 *
 *         EULER EQUATION 1D
 *
 *
 *
 * -----------------------------------------------
 */ 
template <typename Scalar>
MASA::euler_chem_1d<Scalar>::euler_chem_1d()
{
  this->mmsname = "euler_chem_1d";
  this->dimension=1;

  this->register_var("R",&R);
  this->register_var("Cf1_N",&Cf1_N);
  this->register_var("Cf1_N2",&Cf1_N2);

  this->register_var("etaf1_N",&etaf1_N);
  this->register_var("etaf1_N2",&etaf1_N2);

  this->register_var("Ea_N",&Ea_N);
  this->register_var("Ea_N2",&Ea_N2);

  // achtung: need to make into function pointer!
  this->register_var("Function_to_Calculate_K",&Function_to_Calculate_K);

  this->register_var("R_N",&R_N);
  this->register_var("R_N2",&R_N2);

  this->register_var("theta_v_N2",&theta_v_N2);
  this->register_var("M_N",&M_N);

  this->register_var("h0_N",&h0_N);
  this->register_var("h0_N2",&h0_N2);
  this->register_var("K",&K);

}

template <typename Scalar>
int MASA::euler_chem_1d<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("R",1.01);
  err += this->set_var("Cf1_N",1.01);
  err += this->set_var("Cf1_N2",1.01);

  err += this->set_var("etaf1_N",1.01);
  err += this->set_var("etaf1_N2",1.01);

  err += this->set_var("Ea_N",1.01);
  err += this->set_var("Ea_N2",1.01);

  // achtung: need to make into function pointer!
  err += this->set_var("Function_to_Calculate_K",1.01);

  err += this->set_var("R_N",1.01);
  err += this->set_var("R_N2",1.01);

  err += this->set_var("theta_v_N2",1.01);
  err += this->set_var("M_N",1.01);

  err += this->set_var("h0_N",1.01);
  err += this->set_var("h0_N2",1.01);
  err += this->set_var("K",1.01);

  return err;
}

// ----------------------------------------
//   Source Terms
// ----------------------------------------

template <typename Scalar>
Scalar MASA::euler_chem_1d<Scalar>::eval_q_rho_u(Scalar x)
{
  return 0;
}

template <typename Scalar>
Scalar MASA::euler_chem_1d<Scalar>::eval_q_rho_e(Scalar x)
{
  return 0;
}

template <typename Scalar>
Scalar MASA::euler_chem_1d<Scalar>::eval_q_rho_N(Scalar x)
{
  return 0;
}

template <typename Scalar>
Scalar MASA::euler_chem_1d<Scalar>::eval_q_rho_N2(Scalar x)
{
  return 0;
}

// ----------------------------------------
//   Analytical Terms
// ----------------------------------------

template <typename Scalar>
Scalar MASA::euler_chem_1d<Scalar>::eval_exact_t(Scalar x)
{
  return 0;
}

template <typename Scalar>
Scalar MASA::euler_chem_1d<Scalar>::eval_exact_u(Scalar x)
{
  return 0;
}

template <typename Scalar>
Scalar MASA::euler_chem_1d<Scalar>::eval_exact_rho_u(Scalar x)
{
  return 0;
}

template <typename Scalar>
Scalar MASA::euler_chem_1d<Scalar>::eval_exact_rho_N(Scalar x)
{
  return 0;
}

template <typename Scalar>
Scalar MASA::euler_chem_1d<Scalar>::eval_exact_rho_N2(Scalar x)
{
  return 0;
}


// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::euler_chem_1d);

