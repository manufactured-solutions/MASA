// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
// $Author: nick $
// $Id: cns.cpp 13034 2010-09-08 16:17:42Z nick $
//
// rans_sa.cpp: These are the MASA class member functions and constructors
//              For a Reynold Averaged Navier Stokes (RANS) model
//              The Spelart Alamaras
//              For Channel Flow only
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h> 

using namespace MASA;

/* ------------------------------------------------
 *
 *         Spelart Alamaras
 *
 *
 * -----------------------------------------------
 */ 

MASA::rans_sa::rans_sa()
{
  mmsname = "rans_sa";
  dimension=1;

  register_var("Gamma",&Gamma);
  register_var("mu",&mu);

}// done with constructor

void MASA::rans_sa::init_var()
{

  // currently randomly generated -- these are placeholders for now!
  set_var("Gamma",1.01);
  set_var("mu",1.38);

} // done with variable initializer
