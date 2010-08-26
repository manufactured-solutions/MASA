// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// $Author$
// $Id$
//
// euler3d.cpp :program that tests masa against known source term
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <masa.h>
#include <math.h>

using namespace MASA;
using namespace std;

const double pi = acos(-1);
const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

double anQ_p (double x,double y,double z,double p_0,double p_x,double p_y,double p_z,double a_px,double a_py,double a_pz,double L)
{
  double p_an = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L);
  return p_an;
}
 
double anQ_u (double x,double y,double z,double u_0,double u_x,double u_y,double u_z,double a_ux,double a_uy,double a_uz,double L)
{
  double u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L);  
  return u_an;
} 
 
double anQ_v (double x,double y,double z,double v_0,double v_x,double v_y,double v_z,double a_vx,double a_vy,double a_vz,double L)
{
  double v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L);
  return v_an;
}

double anQ_w (double x,double y,double z,double w_0,double w_x,double w_y,double w_z,double a_wx,double a_wy,double a_wz,double L)
{
  double w_an = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L);
  return w_an;
}

double anQ_rho (double x,double y,double z,double rho_0,double rho_x,double rho_y,double rho_z,double a_rhox,double a_rhoy,double a_rhoz,double L)
{ 
  double rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L);
  return rho_an;
}

double SourceQ_e (
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


double SourceQ_u ( // should be 38
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


double SourceQ_v ( // 38
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

double SourceQ_w ( // 38
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

double SourceQ_rho( // 39
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
{
  //variables
  double u_0;
  double u_x;
  double u_y;
  double u_z;
  double v_0;
  double v_x;
  double v_y;
  double v_z;
  double w_0;
  double w_x;
  double w_y;
  double w_z;
  double rho_0;
  double rho_x;
  double rho_y;
  double rho_z;
  double p_0;
  double p_x;
  double p_y;
  double p_z;
  double a_px;
  double a_py;
  double a_pz;
  double a_rhox;
  double a_rhoy;
  double a_rhoz;
  double a_ux;
  double a_uy;
  double a_uz;
  double a_vx;
  double a_vy;
  double a_vz;
  double a_wx;
  double a_wy;
  double a_wz;
  double mu;
  double Gamma;
  double L;    

  // parameters
  double x;
  double y;
  double z;

  // solutions
  double ufield,ufield2,ufield3;
  double vfield,vfield2,vfield3;
  double wfield,wfield2;
  double efield,efield2,efield3;
  double rho,rho2;

  double u_an,u_an2;
  double v_an,v_an2;
  double w_an,w_an2;
  double p_an,p_an2;
  double rho_an,rho_an2;

  // initalize
  int nx = 33;             // number of points
  int ny = 12;  
  int nz = 45;  
  int lx=1;                // length
  int ly=2; 
  int lz=3;   
  double dx=double(lx)/double(nx); // spacing
  double dy=double(ly)/double(ny);
  double dz=double(lz)/double(nz);

  masa_init("euler-test","euler_3d");

  // set params
  masa_init_param();

  // get vars for comparison
  masa_get_param("u_0",&u_0);
  masa_get_param("u_x",&u_x);
  masa_get_param("u_y",&u_y);
  masa_get_param("u_z",&u_z);

  masa_get_param("v_0",&v_0);
  masa_get_param("v_x",&v_x);
  masa_get_param("v_y",&v_y);
  masa_get_param("v_z",&v_z);

  masa_get_param("w_0",&w_0);
  masa_get_param("w_x",&w_x);
  masa_get_param("w_y",&w_y);
  masa_get_param("w_z",&w_z);

  masa_get_param("rho_0",&rho_0);
  masa_get_param("rho_x",&rho_x);
  masa_get_param("rho_y",&rho_y);
  masa_get_param("rho_z",&rho_z);

  masa_get_param("p_0",&p_0);
  masa_get_param("p_x",&p_x);
  masa_get_param("p_y",&p_y);
  masa_get_param("p_z",&p_z);

  masa_get_param("a_px",&a_px);
  masa_get_param("a_py",&a_py);
  masa_get_param("a_pz",&a_pz);

  masa_get_param("a_rhox",&a_rhox);
  masa_get_param("a_rhoy",&a_rhoy);
  masa_get_param("a_rhoz",&a_rhoz);

  masa_get_param("a_ux",&a_ux);
  masa_get_param("a_uy",&a_uy);
  masa_get_param("a_uz",&a_uz);

  masa_get_param("a_vx",&a_vx);
  masa_get_param("a_vy",&a_vy);
  masa_get_param("a_vz",&a_vz);

  masa_get_param("a_wx",&a_wx);
  masa_get_param("a_wy",&a_wy);
  masa_get_param("a_wz",&a_wz);

  masa_get_param("Gamma",&Gamma);
  masa_get_param("mu",&mu);
  masa_get_param("L",&L);

  // check all vars initialized
  masa_sanity_check();

  // evaluate source terms (3D)
  for(int i=0;i<nx;i++)
    for(int j=0;j<ny;j++)    
      for(int k=0;k<nz;k++)
	{
	  x=i*dx;
	  y=j*dy;
	  z=k*dz;

	  //evalulate source terms
	  masa_eval_u_source  (x,y,z,&ufield);
	  masa_eval_v_source  (x,y,z,&vfield);
	  masa_eval_w_source  (x,y,z,&wfield);
	  masa_eval_e_source  (x,y,z,&efield);
	  masa_eval_rho_source(x,y,z,&rho);

	  //evaluate analytical terms
	  masa_eval_u_an        (x,y,z,&u_an);
	  masa_eval_v_an        (x,y,z,&v_an);
	  masa_eval_p_an        (x,y,z,&p_an);
	  masa_eval_rho_an      (x,y,z,&rho_an);	  

	  // check against maple output
	  ufield2   = SourceQ_u  (x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,L);
	  vfield2   = SourceQ_v  (x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,L);
	  wfield2   = SourceQ_w  (x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,L);
	  rho2      = SourceQ_rho(x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,mu,L);
	  efield2   = SourceQ_e  (x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,mu,Gamma,L);
  
	  u_an2     = anQ_u   (x,y,z,u_0,u_x,u_y,u_z,a_ux,a_uy,a_uz,L);
	  v_an2     = anQ_v   (x,y,z,v_0,v_x,v_y,v_z,a_vx,a_vy,a_vz,L);
	  w_an2     = anQ_w   (x,y,z,w_0,w_x,w_y,w_z,a_wx,a_wy,a_wz,L);
	  rho_an2   = anQ_rho (x,y,z,rho_0,rho_x,rho_y,rho_z,a_rhox,a_rhoy,a_rhoz,L);
	  p_an2     = anQ_p   (x,y,z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,L);
	  
	  // test the result is roughly zero
	  ufield3 = fabs(ufield-ufield2);
	  vfield3 = fabs(vfield-vfield2);
	  wfield  = fabs(wfield-wfield2);
	  efield3 = fabs(efield-efield2);
	  rho     = fabs(rho-rho2);

	  u_an   = fabs(u_an-u_an2);
	  v_an   = fabs(v_an-v_an2);
	  rho_an = fabs(rho_an-rho_an2);
	  p_an   = fabs(p_an-p_an2);
	  
	  //masa_display_param();  
	  //cout << endl << ufield << endl << vfield << endl << efield << rho << endl;

	  if(ufield3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "U Field Source Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << ufield3 << endl;
	      cout << "Source term is:                   " << ufield2 << endl;
	      cout << "MASA term is:                     " << ufield << endl;
	      exit(1);
	    }
	  
	  if(u_an > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "U Field Analytical Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << u_an << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(vfield3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "V Field Source Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << vfield3 << endl;
	      cout << "Source term is:                   " << vfield2 << endl;
	      cout << "MASA term is:                     " << vfield << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(v_an > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "V Field Analytical Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << v_an << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(wfield > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "W Field Source Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << wfield << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(w_an > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "W Field Analytical Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << w_an << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(efield3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "Energy Source Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << efield3 << endl;
	      cout << "Source term is:                   " << efield2 << endl;
	      cout << "MASA term is:                     " << efield << endl;
	      cout << x << " " << y << endl;
	      exit(1);
	    }

	  if(p_an > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "P Field Analytical Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << p_an << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(rho > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "RHO Source Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << rho << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(rho_an > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Euler-3d\n";
	      cout << "RHO Analytical Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << rho_an << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }
	  
	}// done iterating

  // tests passed
  return 0;
}
