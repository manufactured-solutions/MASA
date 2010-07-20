//
// program that tests euler-2d against known source term generated from maple
//

#include <masa.h>
#include <math.h>

using namespace MASA;
using namespace std;

const double pi = acos(-1);
const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

double anQ_p (double x,double y,double p_0,double p_x,double p_y,double a_px,double a_py,double L)
{
  double p_an = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L);
  return p_an;
}
  
double anQ_u (double x,double y,double u_0,double u_x,double u_y,double a_ux,double a_uy,double L)
{
  double u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  return u_an;
} 
 
double anQ_v (double x,double y,double v_0,double v_x,double v_y,double a_vx,double a_vy,double L)
{
  double v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);
  return v_an;
}

double anQ_rho (double x,double y,double rho_0,double rho_x,double rho_y,double a_rhox,double a_rhoy,double L)
{ 
  double rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  return rho_an;
}

double SourceQ_e ( // 24
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


double SourceQ_u ( // should be 22
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


double SourceQ_v ( // 22
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


double SourceQ_rho ( // 22
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

{  //variables
  double u_0;
  double u_x;
  double u_y;
  double v_0;
  double v_x;
  double v_y;
  double rho_0;
  double rho_x;
  double rho_y;
  double p_0;
  double p_x;
  double p_y;
  double a_px;
  double a_py;
  double a_rhox;
  double a_rhoy;
  double a_ux;
  double a_uy;
  double a_vx;
  double a_vy;
  double Gamma;
  double mu;
  double L;

  // parameters
  double x;
  double y;

  // solutions -- efield is MASA term, efield2 is maple, efield3 is abs error between them
  double ufield,ufield2,ufield3;
  double vfield,vfield2,vfield3;
  double efield,efield2,efield3;
  double rho,rho2;

  double u_an,u_an2;
  double v_an,v_an2,v_an3;
  double p_an,p_an2;
  double rho_an,rho_an2;

  // initalize
  int nx = 425;  // number of points
  int ny = 768;  
  int lx=3;     // length
  int ly=1; 
  
  double dx=double(lx)/double(nx);
  double dy=double(ly)/double(ny);

  masa_init("euler-test","euler_2d");

  // set params
  masa_init_param();
  
  // get vars
  masa_get_param("u_0",&u_0);
  masa_get_param("u_x",&u_x);
  masa_get_param("u_y",&u_y);
  masa_get_param("v_0",&v_0);
  masa_get_param("v_x",&v_x);
  masa_get_param("v_y",&v_y);

  masa_get_param("rho_0",&rho_0);
  masa_get_param("rho_x",&rho_x);
  masa_get_param("rho_y",&rho_y);

  masa_get_param("p_0",&p_0);
  masa_get_param("p_x",&p_x);
  masa_get_param("p_y",&p_y);

  masa_get_param("a_px",&a_px);
  masa_get_param("a_py",&a_py);

  masa_get_param("a_rhox",&a_rhox);
  masa_get_param("a_rhoy",&a_rhoy);

  masa_get_param("a_ux",&a_ux);
  masa_get_param("a_uy",&a_uy);

  masa_get_param("a_vx",&a_vx);
  masa_get_param("a_vy",&a_vy);

  masa_get_param("Gamma",&Gamma);
  masa_get_param("mu",&mu);
  masa_get_param("L",&L);

  // check that all terms have been initialized
  masa_sanity_check();

  // evaluate source terms (2D)
  for(int i=0;i<nx;i++)
    for(int j=0;j<ny;j++)    
      {
	x=i*dx;
	y=j*dy;
	
	//evalulate source terms
	masa_eval_u_source  (x,y,&ufield);
	masa_eval_v_source  (x,y,&vfield);
	masa_eval_e_source  (x,y,&efield);
	masa_eval_rho_source(x,y,&rho);

	//evaluate analytical terms
	masa_eval_u_an        (x,y,&u_an);
	masa_eval_v_an        (x,y,&v_an);
	masa_eval_p_an        (x,y,&p_an);
	masa_eval_rho_an      (x,y,&rho_an);
	  
	// check against maple
	ufield2 = SourceQ_u   (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,L);
	vfield2 = SourceQ_v   (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,L);
	rho2    = SourceQ_rho (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,L);  
	efield2 = SourceQ_e   (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,Gamma,mu,L);
	
	u_an2   = anQ_u   (x,y,u_0,u_x,u_y,a_ux,a_uy,L);
	v_an2   = anQ_v   (x,y,v_0,v_x,v_y,a_vx,a_vy,L);
	rho_an2 = anQ_rho (x,y,rho_0,rho_x,rho_y,a_rhox,a_rhoy,L);
	p_an2   = anQ_p   (x,y,p_0,p_x,p_y,a_px,a_py,L);
	
	// test the result is roughly zero
	ufield3 = fabs(ufield-ufield2);
	vfield3 = fabs(vfield-vfield2)/fabs(vfield2); // converting to relative error
	efield3 = fabs(efield-efield2)/fabs(efield2); // converting to relative error
	rho    = rho-rho2;
	
	u_an   = u_an-u_an2;
	v_an3  = fabs(v_an-v_an2);
	rho_an = rho_an-rho_an2;
	p_an   = p_an-p_an2;

	if(ufield3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Euler-2d\n";
	    cout << "U Field Source Term\n";
	    cout.precision(16);
	    cout << "Exceeded (relative) Threshold by: " << ufield3 << endl;
	    cout << "Source term is:                   " << ufield2 << endl;
	    cout << "MASA term is:                     " << ufield << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(u_an > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Euler-2d\n";
	    cout << "U Field Analytical Term\n";
	    cout << "Exceeded Threshold by: " << u_an << endl;
	    cout.precision(16);
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(vfield3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Euler-2d\n";
	    cout << "V Field Source Term\n";
	    cout.precision(16);
	    cout << "Exceeded (relative) Threshold by: " << vfield3 << endl;
	    cout << "Source term is:                   " << vfield2 << endl;
	    cout << "MASA term is:                     " << vfield << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(v_an3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Euler-2d\n";
	    cout << "V Field Analytical Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << v_an3 << endl;
	    cout << "Source term is:        " << v_an2 << endl;
	    cout << "MASA term is:          " << v_an << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(efield3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Euler-2d\n";
	    cout << "Energy Source Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << efield3 << endl;
	    cout << "Source term is:        " << efield2 << endl;
	    cout << "MASA term is:          " << efield << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(p_an > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Euler-2d\n";
	    cout << "P Field Analytical Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << p_an << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(rho > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Euler-2d\n";
	    cout.precision(16);
	    cout << "RHO Source Term\n";
	    cout << "Exceeded Threshold by: " << rho << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(rho_an > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Euler-2d\n";
	    cout.precision(16);
	    cout << "RHO Analytical Term\n";
	    cout << "Exceeded Threshold by: " << rho_an << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

      } // done iterating
  // tests passed

}
