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
// $Author$
// $Id$
//
// c_euler_example.c: this is an example of the heat equation in 1d for C
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa.h>
#include <stdio.h>

int main()
{
  double sol;
  
  // init
  cmasa_init("nick","euler_1d");
  cmasa_init_param();

  // list
  cmasa_list_mms();
  cmasa_display_param();

  //check all initialized properly
  cmasa_sanity_check();
  sol=cmasa_eval_1d_source_u(1.2);
  printf("\nt source: %g\n",sol);

  return 0; // done
}
