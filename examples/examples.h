// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011 The PECOS Development Team
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
// $Id: euler_example.cpp 19214 2011-03-29 19:33:16Z karl $
//
// euler_example.cpp:
// this is an example of the API used for calling the 2D euler equation
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef __cplusplus
#include <limits>
#include <iostream>

template <typename Scalar>
void test(Scalar input)
{
  Scalar MASA_VAR_DEFAULT = -12345.67;
  Scalar uninit = -1.33;
  Scalar thresh = 5 * std::numeric_limits<Scalar>::epsilon();

  if( fabs((input - MASA_VAR_DEFAULT)/MASA_VAR_DEFAULT) < thresh)
    {
      exit(1);
    }

  if(fabs((input - uninit)/uninit) < thresh)
    {
      exit(1);
    }
  
}

template <typename Scalar>
Scalar temp_function(Scalar T);

template<typename Scalar>
Scalar funct(Scalar T);

#endif // __cplusplus

#ifdef __cplusplus
extern "C" {
#endif

void test(double input)
{
  double MASA_VAR_DEFAULT = -12345.67;
  double uninit = -1.33;
  double thresh = 5 * 1e-15;

  if( fabs((input - MASA_VAR_DEFAULT)/MASA_VAR_DEFAULT) < thresh)
    {
      exit(1);
    }

  if(fabs((input - uninit)/uninit) < thresh)
    {
      exit(1);
    }
  
}

#ifdef __cplusplus
}
#endif
