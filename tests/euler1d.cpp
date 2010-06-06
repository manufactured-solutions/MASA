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
  // parameters
  double param = 1.2;
  double x=.5;

  // solutions
  double ufield,ufield2;
  double efield,efield2;
  double rho,rho2;

  // initalize
  masa_init("euler-test","euler_1d");

  // set params
  masa_set_param("u_0",param);
  double u_0=param;
  masa_set_param("u_x",param);
  double u_x=param;

  masa_set_param("rho_0",param);
  double rho_0=param;
  masa_set_param("rho_x",param);
  double rho_x=param;

  masa_set_param("p_0",param);
  double p_0=param;
  masa_set_param("p_x",param);
  double p_x=param;

  masa_set_param("a_px",param);
  double a_px=param;

  masa_set_param("a_rhox",param);
  double a_rhox=param;

  masa_set_param("a_ux",param);
  double a_ux=param;

  masa_set_param("Gamma",param);
  double Gamma=param;
  masa_set_param("mu",param);
  double mu=param;
  masa_set_param("L",param);
  double L=param;

  // evaluate source terms (1D)
  masa_sanity_check();
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
