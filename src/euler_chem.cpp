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
}

template <typename Scalar>
int MASA::euler_chem_1d<Scalar>::init_var()
{
  int err = 0;
  err += this->set_var("R",1.01);

  return err;
}

// ----------------------------------------
//   Source Terms
// ----------------------------------------



// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::euler_chem_1d);

