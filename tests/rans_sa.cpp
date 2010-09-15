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
// $Id: euler1d.cpp 13404 2010-09-15 02:56:33Z nick $
//
// rans_sa.cpp :program that tests masa against known source term
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

//
// OBVIOUSLY NOT FINISHED!!
//

#include <config.h>
#include <masa.h>
#include <math.h>

#include <iostream>
#include <stdlib.h>

using namespace MASA;
using namespace std;

const double pi = acos(-1);
const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

int main()
{
  // initalize
  masa_init("spelart alamaras test","rans_sa");

  // initialize the default parameters
  masa_init_param();

  // check that all terms have been initialized
  masa_sanity_check();


  return 0;
}

