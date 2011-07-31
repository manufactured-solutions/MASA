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
// $Id$
//--------------------------------------------------------------------------

/*!
 * \file laplacian.c
 *
 * \brief A simple standalone tutorial program which solves a
 * Laplacian using finite-differening and uses the MASA library to aid
 * in verification.  For clarity, this program is self contained with
 * all required routines present in two source files (and one header
 * file). Note that for convenience, we use naive data structures and
 * form a full dense matrix. In addition, no particular performance
 * considerations are applied.
 *
 */

#include "laplacian.h"

int main(int argc, char *argv[])
{
  int n;
  double length;
  pstruct model;		/* primary model data structure */

  /* Parse command-line */

  if(argc < 2)
    {
      printf("\nUsage: laplacian [num_pts] [length]\n\n");
      printf("where \"num_pts\" is the desired number of mesh points and \n");
      printf("\"length\" is the physical length-scale dimension in one direction\n\n");
      exit(1);
    }
  else
    {
      n      = atoi(argv[1]);
      length = (double) atof(argv[2]);
    }
  
  /* Problem Initialization */

  problem_initialize (n,length,&model);  
  assemble_matrix    (1,&model);		
  init_masa          (&model);
  apply_bcs          (&model);

  /* Solve */

  solve_gauss       (&model);

  /* Compute Error */

  printf("\n** Error Analysis\n");
  printf("   --> npts     = %i\n",model.npts);
  printf("   --> h        = %12.5e\n",model.h);
  printf("   --> l2 error = %12.5e\n",compute_l2_error(&model));

  return 0;
}

