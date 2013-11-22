//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012,2013 The PECOS Development Team
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

#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<masa.h>
#include<math.h>

/* Primary model data-structure */

typedef struct pstruct {
  double  *phi;			/*!< solution variable                    */
  double  *rhs;			/*!< right-hand side forcing function     */
  double **A;			/*!< linear system matrix                 */
  double   h;			/*!< mesh sizing                          */
     int   n;			/*!< problem size                         */
     int   npts;		/*!< number of points in single direction */
     int   pad;			/*!< pad dimension for ghost points       */
} pstruct;

/* Supported finite-difference stencils */

enum fd_types
  {
    central_2nd_order,
    central_4th_order
  };

/* Function prototypes */

  void apply_bcs            (pstruct *model);
  void assemble_matrix      (const int fd_method, pstruct *model);
  void compute_error        (pstruct *model);
double compute_l2_error     (pstruct *model);
   int converged            (double *a, double *b, double eps, const int n, double *diff);
  void enforce_dirichlet_bc (const int col_id, const int row_id, const double value, pstruct *model);
  void init_masa            (pstruct *model);
  void print_matrix         (pstruct *model);
  void problem_initialize   (const int, const double length, pstruct *model);
  void solve_gauss          (pstruct *model);
  void solve_cg             (pstruct *model);

/* conjugate gradient */

int cg(int n,double* A,double* b,double* x,double* residual,int* iterations);
