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
// $Id: c_euler1d.c 13404 2010-09-15 02:56:33Z nick $
//
// sod.cpp : test sod mms
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using namespace MASA;
using namespace std;

int main()
{
  double out;
  double x=1.1;
  double t=2.3;

  masa_init("sod-test","sod_1d");
  masa_init_param();
  masa_sanity_check();
  
  out = masa_eval_rho_source(x,t);
  out = masa_eval_rho_u_source(x,t);

}
