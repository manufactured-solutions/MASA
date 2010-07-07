
#include <iostream>
#include <math.h>

using namespace std;

// takes two vectors, returns scalar product
double scalar(double *x1,double *x2, int n)
{
  double scal=0;
  for(int i=0;i<n;i++)
    scal += x1[i]*x2[i];
  
  return scal;
}

// return L infinity norm
double linf(double *x1,int n)
{
  double scal=0;
  for(int i=0;i<n;i++)
    {    
      //cout << "hit " << x1[i] << endl;
	  
      if(x1[i]>scal)
	{
	  scal=x1[i];
	}
    }
  return scal;
}

// return L2 norm
double l2(double *x1,int n)
{
  double scal=0;
  for(int i=0;i<n;i++)
    {    
      scal += x1[i]*x1[i];
    }
  //scal=sqrt(scal);
  return scal;
}

// j is column
// i is row

// arr[i][j] = arr[i*ncolumn + j]
//  x[i]     = x  [i]
int cg(int n,double* A,double* b,double* x)
{
  double r [n]; // residual
  double ro[n]; // old residual
  double p[n];  // search vector
  double alpha; // step length
  double beta;  // improvement
  
  double thresh = 1e-15;
  int it = 0;
  int count = 20;
  double temp[n];  // temp vector
  
  for(int i=0;i<n;i++)
    {
      // initialize residual to RHS
      r [i] = b[i];
      ro[i] = b[i];
      // init search vector to RHS
      p[i] = b[i];

      // init solution vector to zero
      x[i] = 0;

      temp[i] = 0;
    }

  //begin iteration
  //while(linf(&r[0],n)>thresh)
  while(n > it)
    {
      // calc A p_{n-1}
      for(int ii=0;ii<n;ii++)
	{
	  temp[ii]=0;
	  for(int jj=0;jj<n;jj++)
	    {
	      temp[ii] += A[ii*n+jj]*p[jj];
	    }
	}
      
      // update step length
      alpha = scalar(&ro[0],&ro[0],n)/scalar(&p[0],&temp[0],n);
      //cout << "alpha is: " << alpha << endl;

      // update approx solution & residual 
      for(int ii=0;ii<n;ii++)
	{
	  x[ii] +=  alpha*p[ii];	  
	  r[ii] = ro[ii]-alpha*temp[ii];
	}
	
      // calc improvement this step
      beta = scalar(&r[0],&r[0],n)/scalar(&ro[0],&ro[0],n);

      // update search direction
      for(int ii=0;ii<n;ii++)
	{
	  p[ii]  = r[ii] + beta * p[ii];
	  ro[ii] = r[ii]; // also update old residual
	}
      
      it++;
      //cout << r[0] << endl;
      //cout << "current thresh exceeded by " << l2(&r[0],n) << endl;

    }//done with iteration

  cout << "Performed " << it << " Iteration(s)" << endl;
  
}

// this is just a test routine to check the conjugate gradient solver
//

int main()
{
  int n=3;        // size of everything
  double A[n][n]; // matrix
  double b[n];    // RHS
  double x[n];    // solution vector
   
  // init n:
  for(int i=0;i<n;i++)
    {
      /*
      for(int j=0;j<n;j++)
	{
	  A[i][j]=i*n+j+1;	  
	}
      */

      A[0][0]=1;
      A[0][1]=2;
      A[0][2]=3;

      A[1][0]=2;
      A[1][1]=1;
      A[1][2]=5;

      A[2][0]=3;
      A[2][1]=5;
      A[2][2]=1;

      b[0]=1;
      b[1]=2;
      b[2]=3;
      //x[i]=0;
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

  cout << "\n\nb is: \n";
  for(int i=0;i<n;i++)
    {
      cout << "[ " << b[i] << " ]\n";
    }

  cout << "\n\nx is: \n";
  for(int i=0;i<n;i++)
    {
      cout << "[ " << x[i] << " ]\n";
    }

}// end program
