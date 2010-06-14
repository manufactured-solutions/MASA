//
// program that tests euler-1d against known source term generated from maple
//

#include <masa.h>
#include <math.h>

using namespace MASA;
using namespace std;

const double pi = acos(-1);
const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

double anQ_p (double x,double p_0,double p_x,double a_px,double L)
{
  double p_an = p_0 + p_x * cos(a_px * pi * x / L);
  return p_an;
}
  
double anQ_u (double x,double u_0,double u_x,double a_ux,double L)
{
  double u_an = u_0 + u_x * sin(a_ux * pi * x / L);
  return u_an;
} 
 
double anQ_rho (double x,double rho_0,double rho_x,double a_rhox,double L)
{ 
  double rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  return rho_an;
}

double SourceQ_e ( // 12
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double);

double SourceQ_u ( // should be 10
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double);

double SourceQ_rho ( // 10
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double);

int main()
{
  //variables 
  double u_0;
  double u_x;
  double rho_0;
  double rho_x;
  double p_0;
  double p_x;
  double a_px;
  double a_rhox;
  double a_ux;
  double Gamma;
  double mu;
  double L;

  // parameters
  double x;

  //problem size
  int nx = 2000;  // number of points
  int lx=10;     // length
  double dx=double(lx)/double(nx);

  // solutions
  double ufield,ufield2;
  double efield,efield2;
  double rho,rho2;

  double u_an,u_an2;
  double p_an,p_an2;
  double rho_an,rho_an2;

  // initalize
  masa_init("euler-test","euler_1d");

  // initialize the default parameters
  masa_init_param();

  // get defaults for comparison to source terms
  masa_get_param("u_0",&u_0);
  masa_get_param("u_x",&u_x);

  masa_get_param("rho_0",&rho_0);
  masa_get_param("rho_x",&rho_x);

  masa_get_param("p_0",&p_0);
  masa_get_param("p_x",&p_x);

  masa_get_param("a_px",&a_px);
  masa_get_param("a_rhox",&a_rhox);
  masa_get_param("a_ux",&a_ux);

  masa_get_param("Gamma",&Gamma);
  masa_get_param("mu",&mu);
  masa_get_param("L",&L);

  // check that all terms have been initialized
  masa_sanity_check();

  // evaluate source terms (1D)
  for(int i=0;i<nx;i++)
    {
      x=i*dx;
      
      // ask for masa answer
      masa_eval_u_source  (x,&ufield);
      masa_eval_e_source  (x,&efield);
      masa_eval_rho_source(x,&rho);

      //evaluate analytical terms
      masa_eval_u_an        (x,&u_an);
      masa_eval_p_an        (x,&p_an);
      masa_eval_rho_an      (x,&rho_an);

      // get fundamental source term solution
      ufield2   = SourceQ_u  (x,u_0,u_x,rho_0,rho_x,p_0,p_x,a_px,a_rhox,a_ux,L);
      rho2      = SourceQ_rho(x,u_0,u_x,rho_0,rho_x,p_0,p_x,a_px,a_rhox,a_ux,L);
      efield2   = SourceQ_e  (x,u_0,u_x,rho_0,rho_x,p_0,p_x,a_px,a_rhox,a_ux,Gamma,mu,L);
  
      u_an2   = anQ_u   (x,u_0,u_x,a_ux,L);
      rho_an2 = anQ_rho (x,rho_0,rho_x,a_rhox,L);
      p_an2   = anQ_p   (x,p_0,p_x,a_px,L);

      // test the result is roughly zero
      ufield=ufield-ufield2;
      efield=efield-efield2;
      rho   =rho-rho2;

      u_an   = u_an-u_an2;
      rho_an = rho_an-rho_an2;
      p_an   = p_an-p_an2;
   
      if(ufield > threshold)
	{
	  cout << "\nMASA REGRESSION TEST FAILED: Euler-1d\n";
	  cout << "U Field Source Term\n";
	  cout << "Exceeded Threshold by: " << ufield << endl;
	  exit(1);
	}

      if(u_an > threshold)
	{
	  cout << "\nMASA REGRESSION TEST FAILED: Euler-1d\n";
	  cout << "U Field Analytical Term\n";
	  cout << "Exceeded Threshold by: " << u_an << endl;
	  exit(1);
	}

      if(efield > threshold)
	{
	  cout << "\nMASA REGRESSION TEST FAILED: Euler-1d\n";
	  cout << "Energy Source Term\n";
	  cout << "Exceeded Threshold by: " << efield << endl;
	  exit(1);
	}

      if(p_an > threshold)
	{
	  cout << "\nMASA REGRESSION TEST FAILED: Euler-1d\n";
	  cout << "P Field Analytical Term\n";
	  cout << "Exceeded Threshold by: " << p_an << endl;
	    exit(1);
	}
      
      if(rho > threshold)
	{
	  cout << "\nMASA REGRESSION TEST FAILED: Euler-1d\n";
	  cout << "RHO Source Term\n";
	  cout << "Exceeded Threshold by: " << rho << endl;
	  exit(1);
	}
      
      if(rho_an > threshold)
	{
	  cout << "\nMASA REGRESSION TEST FAILED: Euler-1d\n";
	  cout << "RHO Analytical Term\n";
	  cout << "Exceeded Threshold by: " << rho_an << endl;
	  exit(1);
	}

    } // done interating 
  // tests passed

}
