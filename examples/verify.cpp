// this is just a test routine to check the conjugate gradient solver
//
#include <iostream>
#include <math.h>
#include <masa.h>

using namespace MASA;
using namespace std;

double a = 3.45;
double b = -1.89;
//double c = .071;
//double d = 12.407;

double c = 0;
double d = 0;

// used to solve the system 
int cg(int,double*,double*,double*);

/*
 *
 * This routine takes the spacial location as input,
 * and returns the value of the polynomial (after differentiation) 
 * at that location.
 *
 */
double source(double x)
{
  double ret_val = (20 * a * x*x*x) + (6 * b * x);
  return ret_val;
}

/*
 *
 * This routine takes the spacial location as input,
 * and returns the value of the polynomial (solution)
 * at that location.
 *
 */
double sol(double x)
{
  double ret_val = (a* x*x*x*x*x) + (b* x*x*x) + (c*x) + d;
  return ret_val;
}

/*
 *
 * This routine takes the number of points as input,
 * and returns the value of a 2nd order finite difference
 * discretization for the 2nd derivative operator
 * 
 */
int finite_diff(int n,double L)
{
  // n is the size of the system (in number of points)
  // L is the length of the problem (starting from zero)

  double A[n][n]; // matrix
  double b[n];    // RHS
  double x[n];    // solution vector

  double residual=0;
   
  double dx = double(L/n); // this is also known as 'h'
  int lto;

  // init n:
  for(int i=0;i<n;i++)
    {
      for(int j=0;j<n;j++)
	{
	  lto = j-i;
	  switch(lto)
	    {
	    case (-1): // one less than diagonal
	      A[i][j] = 1;
	      break;
	      
	    case (0): // diagonal
	      A[i][j] = -2;
	      break;
	      
	    case (1): //  one more than diagonal
	      A[i][j] = 1;
	      break;

	    default:    // not near diagonal
	      A[i][j]=0; // set to zero
	      break;
	      
	    } // done with case construct	  
	}
    } // done building matrix A

  
  for(int i=0;i<n;i++) // build solution matrix
    {      
      b[i]=source(i*dx); // this is the analytical solution
    }
  
  // assume (n x n) (n x 1) and (n x 1)
  cg(n,&A[0][0],&b[0],&x[0]);
  
  cout << "A is: \n";
  for(int i=0;i<n;i++)
    {
      for(int j=0;j<n;j++)
	{
	  // cout << "[ " << i << " " << j <<  " " << x[i][j] << " ] ";
	  cout << "[ " << A[i][j] << " ] ";
	}
      cout << endl;
    }

  cout << "\n\nb (RHS) is: \n";
  for(int i=0;i<n;i++)
    {
      cout << i*dx <<  " " << b[i] << "\n";
    }
  
  cout << "\n\nx (solution) is: \n";
  for(int i=0;i<n;i++)
    {
      cout << i*dx << " " << x[i] << "\n";
    }

  cout << "\n\nanalytical solution is: \n";
  for(int i=0;i<n;i++)
    {
      cout << i*dx << " " << sol(i*dx) << endl;
    }


  // now lets check this is in fact the correct solution
  for(int i;i<n;i++)
    {
      residual += (dx*dx*x[i])-sol(i*dx);
    }

  cout << "residual is: " << residual << endl;

}// end program

int main()
{

  // Here:
  // We calculate the solution using a 2nd order finite difference
  // Then:
  // We compare this to the analytical result stored in the masa library
  // And observe the cnvergence behaviour as we refine the mesh.

  // start up masa
  double tempx,tfield;

  /*  masa_init("polynomial example","poly_example");
  masa_init_param();
  masa_sanity_check();

  finite_diff(10,2);
  masa_eval_t_source  (tempx,&tfield);
  
  finite_diff(20,2);

  finite_diff(40,2);
  */

  finite_diff(10,.5);

}
