
/*!
 * \file laplacian.c
 *
 * \brief A simple standalone program which solves a Laplacian using
 * finite-differening and uses the MASA library to aid in
 * verification.  For clarity, this program is self contained with all
 * required routines present in a single source file. Note that for
 * convenience, we use naive data structures, which do not account for
 * the sparse matrix structure and no particular performance
 * considerations are applied.
 * 
 * \author Karl W. Schulz (karl@ices.utexas.edu)
 *
 */

#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<masa.h>
#include<math.h>

typedef struct pstruct {
  double  *phi;			/*!< solution variable                    */
  double  *rhs;			/*!< right-hand side forcing function     */
  double **A;			/*!< linear system matrix                 */
  double   h;			/*!< mesh sizing                          */
  int      n;			/*!< problem size                         */
  int   npts;			/*!< number of points in single direction */
  int    pad;			/*!< pad dimension for ghost points       */
} pstruct;

void problem_initialize (const int, const double length, pstruct *model);
void apply_bcs          (pstruct *model);
void assemble_matrix    (int fd_method, pstruct *model);
void init_masa          (pstruct *model);
void compute_error      (pstruct *model);
void print_matrix       (pstruct *model);
int  converged          (double *a, double *b, double eps, int n, double *diff);
void apply_dirichlet_bc (const int col_id, const int row_id, const double value, pstruct *model);
void solve_gauss        (pstruct *model);
double compute_l2_error (pstruct *model);

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

void compute_error(pstruct *model)
{
  double error = 0.0;
  double xval, yval;
  int i,j;
  int index;

  for(i=0;i<model->npts;i++)
    {
      for(j=0;j<model->npts;j++)
	{
	  index = j+(i*model->npts);
	  
	  xval = (i)*model->h;
	  yval = (j)*model->h;
	  
	  printf("%3i: Num solution = %12.5e, Analytic solution = %12.5e (diff = %f)\n",index,
		 model->phi[index], masa_eval_2d_exact_phi(xval,yval),
		 fabs(model->phi[index]-masa_eval_2d_exact_phi(xval,yval)));
	}
      printf("\n");
    }
  return;
}

void init_masa(pstruct *model)
{

  int i,j;
  int index;
  double xval, yval;

  printf("\n** Initializing MASA\n");

  masa_init("A C Laplacian Example","laplace_2d");
  
  printf("   --> applying RHS forcing function\n");

  /* Define RHS for solution domain */

  for(i=0;i<model->npts;i++)
    {
      for(j=0;j<model->npts;j++)
	{
	  xval  = (i)*model->h;
	  yval  = (j)*model->h;
	  index = j+(i*model->npts);

	  model->rhs[index] = masa_eval_2d_source_f(xval,yval);
	}
    }
  printf("   --> complete\n");
  return;
}

void print_matrix(pstruct *model)
{
  int i,j;

  printf("     j =");
  for(j=0;j<model->n;j++)
    printf("%6i ",j);
  printf("\n\n");
  
  for(i=0;i<model->n;i++)
    {
      printf("i = %2i: ",i);
      for(j=0;j<model->n;j++)
	printf("%6.2f ",model->A[i][j]);
      printf("\n");
    }

  printf("\nRHS:\n");
  for(j=0;j<model->n;j++)
    printf("%6.2f\n",model->rhs[j]);

}

void apply_dirichlet_bc(const int row_id, const int col_id, const double value, pstruct *model)
{
  int i,j;
  int index;

  index = col_id+(row_id*model->npts);

  for(j=0;j<model->n;j++)
    model->A[index][j] = 0.0;

  model->A[index][index] = 1.0;
  model->rhs[index]      = value;

  return;
}

void assemble_matrix(int fd_method, pstruct *model)
{
  int i,j;
  int index,index2;
  const double h_squared = model->h*model->h;

  const int n = model->n;

  printf("\n** Assembling linear system\n");

  for(i=0;i<model->npts;i++)
    for(j=0;j<model->npts;j++)
      {
	index = j+(i*model->npts);

	model->A[index][index] = -4.0/h_squared;

	if(index > 0)
	  model->A[index-1][index] = 1.0/h_squared;

	if(index < model->n-1 )
	  model->A[index+1][index] = 1.0/h_squared;

	if(i > 0)
	  {
	    index2 = j+(i-1)*model->npts;
	    model->A[index2][index] = 1.0/h_squared;
	  }

	if(i < model->npts-1)
	  {
	    index2 = j+(i+1)*model->npts;
	    model->A[index2][index] = 1.0/h_squared;
	  }
      }

  printf("   --> assembly complete\n");
  return;
}

void apply_bcs(pstruct *model)
{
  int i,j;
  int index;
  double soln;
  const int n    = model->n;
  const int npts = model->npts;

  double xval,yval;

  assert(model->pad == 1);

  if(model->pad >= 1)
    {
      i=0;                              /* BCs for north boundary */
      for(j=0;j<model->npts;j++)      
	{
	  xval = (i)*model->h;
	  yval = (j)*model->h;
	  soln = masa_eval_2d_exact_phi(xval,yval);
	  
	  apply_dirichlet_bc(i,j,soln,model);
	}

      i=model->npts-1;                 /* BCs for south boundary */
      for(j=0;j<model->npts;j++)      
	{
	  xval = (i)*model->h;
	  yval = (j)*model->h;
	  soln = masa_eval_2d_exact_phi(xval,yval);

	  apply_dirichlet_bc(i,j,soln,model);
	}

      j=0;                            /* BCs for west boundary */
      for(i=0;i<model->npts;i++)      
	{
	  xval = (i)*model->h;
	  yval = (j)*model->h;
	  soln = masa_eval_2d_exact_phi(xval,yval);

	  apply_dirichlet_bc(i,j,soln,model);
	}

      j=model->npts-1;                /* BCs for east boundary */
      for(i=0;i<model->npts;i++)      
	{
	  xval = (i)*model->h;
	  yval = (j)*model->h;
	  soln = masa_eval_2d_exact_phi(xval,yval);

	  apply_dirichlet_bc(i,j,soln,model);
	}
    }

  return;
}

void problem_initialize(const int npts, const double length, pstruct *model)
{

/*! \fn           problem_initialize
 *  \brief        Allocates memory for desired problem size.
 *  \param npts    Desired # points in one direction (resulting solution vector is npts*npts)
 *  \param length Dimension of unit square
 *  \param model  Laplacian data structure.
 */

  int i;

  assert(npts   > 1 );
  assert(length > 0.);

  printf("\n** Initializing problem\n");

  model->npts   = npts;
  model->h      = length/(npts-1.0);
  model->n      = npts*npts;
  model->pad    = 1;

  printf("   --> mesh size           = %-12.3f\n",model->h);
  printf("   --> ghost node pad size = %i\n",model->pad);

  /* Perform memory allocation */

  model->rhs = (double *) calloc(sizeof(double  ),     model->n);
  model->phi = (double *) calloc(sizeof(double  ),     model->n);
  model->A   = (double **)calloc(sizeof(double *),     model->n);

  if(model->rhs == NULL || model->phi == NULL || model->A  == NULL )
    {
      printf("ERROR: Unable to allocate memory\n");
      exit(0);
    }

  for(i=0;i<model->n;i++)
    {
      model->A[i] = calloc(sizeof(double),model->n);
      if(model->A[i] == NULL)
	{
	  printf("ERROR: Unable to allocate memory for matrix A\n");
	  exit(0);
	}
    }    
  
  printf("   --> memory initialization complete\n");
  return;
}


void solve_gauss(pstruct *model)
{
  int it=0;
  int i,j, itmp;
  double *phi_old;
  int   **sparse_index;
  int    *sparse_count;
  int col_index;
  double diff            = 0.0;
			
  const double  ZERO_TOL = 1.0e-20;
  const double CONVG_TOL = 1.0e-9;
  const int    ITER_MAX  = 30000;
  const int           n  = model->n;

  printf("\n** Solving system using Gauss-Seidel\n");

  phi_old      =        calloc(sizeof(double),model->n);
  sparse_index = (int **)calloc(sizeof(int *),model->n);
  sparse_count = ( int *)calloc(sizeof(int)  ,model->n);

  assert(phi_old != NULL && sparse_index != NULL && sparse_count != NULL);

  printf("   --> Building sparse matrix map\n");

  /* Build sparse-matrix mapping */

  for(i=0;i<n;i++)
    {

      /* Count number of non-zero entries on this row */

      sparse_count[i] = 0;
      for(j=0;j<n;j++)
	if(fabs(model->A[i][j]) > ZERO_TOL)
	  sparse_count[i]++;

      /* Allocate space and store non-zero entries for this row */

      sparse_index[i] = (int *)calloc(sizeof(int), sparse_count[i]);

      sparse_count[i] = 0;
      for(j=0;j<n;j++)
	if(fabs(model->A[i][j]) > ZERO_TOL)
	  sparse_index[i][sparse_count[i]++] = j;
    }

  /* iterate until convergencd */

  printf("   --> Solving linear system\n");

  while(it < ITER_MAX)
    {
      /* Save last iterate */

      for(i=0;i<model->n;i++)
	phi_old[i] = model->phi[i];

      for(i=0;i<n;i++)		/* loop-over rows  */
	{
	  model->phi[i] = model->rhs[i];
	  
	  for(j=0;j<sparse_count[i];j++)
	    {
	      col_index      = sparse_index[i][j];
	      if(col_index != i)
		model->phi[i] -= model->A[i][col_index]*model->phi[col_index];
	    }
	  
	  model->phi[i] =  model->phi[i] / model->A[i][i];
	}

      it++;

      if(converged(model->phi,phi_old,CONVG_TOL,model->n,&diff)) break;

    } /* end main iterative loop */

  free(phi_old);
  for(i=0;i<model->n;i++)
    free(sparse_index[i]);
  free(sparse_index);
  
  free(sparse_count);

  printf("   --> Converged in %i iters: diff: %15.7g\n", it, diff);
}

int converged(double *a, double *b, double eps, int n, double *diff)
{
  int i;
  double sum=0.;
  for (i=0; i<n; i++)
    sum += (a[i]-b[i])*(a[i]-b[i]);

  sum /= ((double)(n));

  *diff=sqrt(sum);

  if(*diff < eps)
    return(1);
  else
    return(0);
}

double compute_l2_error(pstruct *model)
{
  double l2_error = 0.0;
  double xval,yval;
  double diff;
  int i,j,index;

  for(i=0;i<model->npts;i++)
    for(j=0;j<model->npts;j++)
      {
	index = j+(i*model->npts);
	
	xval = (i)*model->h;
	yval = (j)*model->h;

	diff = masa_eval_2d_exact_phi(xval,yval)-model->phi[index];

	l2_error += diff*diff;
      }

  l2_error = sqrt(l2_error / model->n );

  return(l2_error);
}
