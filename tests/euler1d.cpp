//
// program that tests euler-1d against known source term generated from maple
//

#include <masa.h>

using namespace MASA;
using namespace std;

#define PI = acos(-1)

const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

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
  double x=.5;

  //problem size
  double lx,ly;
  double dx,dy;
  int nx,ny;

  // initialize
  nx = 10;  // number of points
  ny = 10;  
  lx=1;     // length
  ly=1; 

  dx=double(lx/nx);
  dy=double(ly/ny);

  // solutions
  double ufield,ufield2;
  double efield,efield2;
  double rho,rho2;

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
  
  masa_eval_u_source  (x,&ufield);
  masa_eval_e_source  (x,&efield);
  masa_eval_rho_source(x,&rho);

  ufield2   = SourceQ_u  (x,u_0,u_x,rho_0,rho_x,p_0,p_x,a_px,a_rhox,a_ux,L);
  rho2      = SourceQ_rho(x,u_0,u_x,rho_0,rho_x,p_0,p_x,a_px,a_rhox,a_ux,L);
  efield2   = SourceQ_e  (x,u_0,u_x,rho_0,rho_x,p_0,p_x,a_px,a_rhox,a_ux,Gamma,mu,L);
  
  // test the result is roughly zero
  ufield=ufield-ufield2;
  efield=efield-efield2;
  rho   =rho-rho2;
   
  if(ufield > threshold)
    {
      cout << "\nMASA REGRESSION TEST FAILED: Euler-1d\n";
      cout << "U Field Source Term\n";
      cout << "Exceeded Threshold by: " << ufield << endl;
      exit(1);
    }

  if(efield > threshold)
    {
      cout << "\nMASA REGRESSION TEST FAILED: Euler-1d\n";
      cout << "Energy Source Term\n";
      cout << "Exceeded Threshold by: " << efield << endl;
      exit(1);
    }

  if(rho > threshold)
    {
      cout << "\nMASA REGRESSION TEST FAILED: Euler-1d\n";
      cout << "RHO Source Term\n";
      cout << "Exceeded Threshold by: " << rho << endl;
      exit(1);
    }

  // tests passed

}
