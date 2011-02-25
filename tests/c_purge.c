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
// $Author: karl $
// $Id: c_misc.c 13969 2010-09-23 13:55:56Z karl $
//
// c_purge.c :program that tests masa purge function
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double MASA_DEFAULT = -12345.67;

int main()
{
  double u_0;

  // start problem
  masa_init("euler-test","euler_1d");

  u_0 = masa_get_param("u_0");
  if(u_0 == MASA_DEFAULT)
    {
      printf("\nMASA ERROR:: Variables not being auto initalized!\n");
      return 1;
    }
  
  // now purge default values
  masa_purge_default_param();
  // values should be set to default: checking
  u_0 = masa_get_param("u_0");
  if(u_0 != MASA_DEFAULT)
    {
      printf("\nMASA ERROR:: Variables not being purged properly!\n");
      return 1;
    }
  
  // steady as she goes
  return 0;
}
