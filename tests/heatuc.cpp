#include <math.h>
#include <masa.h>

using namespace std;
using namespace MASA;

const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

double SourceQ_t_1d(
  double x,
  double t,
  double A_x,
  double A_t,
  double D_t,
  double k_0,
  double cp_0,
  double rho)
{
  double Q_T = cos(A_x * x + A_t * t) * cos(D_t * t) * k_0 * A_x * A_x - (sin(A_x * x + A_t * t) * cos(D_t * t) * A_t + cos(A_x * x + A_t * t) * sin(D_t * t) * D_t) * rho * cp_0;
  return Q_T;
}

double SourceQ_t_2d (
  double x,
  double y,
  double t,
  double A_x,
  double A_t,
  double B_y,
  double B_t,
  double D_t,
  double rho,
  double k_0,
  double cp_0)
{
  double Q_T = -(sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(D_t * t) * A_t + cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(D_t * t) * B_t + cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(D_t * t) * D_t) * rho * cp_0 + (A_x * A_x + B_y * B_y) * cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(D_t * t) * k_0;
  return Q_T;
}

double SourceQ_t_3d (
  double x,
  double y,
  double z,
  double t,
  double A_x,
  double A_t,
  double B_y,
  double B_t,
  double C_z,
  double C_t,
  double D_t,
  double k_0,
  double cp_0,
  double rho)
{
  double Q_T = -(sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * A_t + cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * B_t + cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(C_z * z + C_t * t) * cos(D_t * t) * C_t + cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * sin(D_t * t) * D_t) * rho * cp_0 + (A_x * A_x + B_y * B_y + C_z * C_z) * cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * k_0;
  return Q_T;
}

int main()
{
  double tfield,tfield2;
  double param=1.2;
  double t= 1;
  double x=.5;
  double y=.4;
  double z=.3;

  /// -----------------------------------------------------------------------
  // initalize 1D
  /// -----------------------------------------------------------------------
  masa_init("temp-test-1d","heateq_1d_unsteady_const");

  masa_set_param("A_x",param);
  double A_x=param;

  masa_set_param("A_t",param);
  double A_t=param;

  masa_set_param("D_t",param);
  double D_t=param;

  masa_set_param("cp_0",param);
  double cp_0=param;

  masa_set_param("rho",param);
  double rho=param;

  masa_set_param("k_0",param);
  double k_0=param;

  // evaluate source terms (1D)
  masa_sanity_check();
  masa_eval_t_source(x,t,&tfield);
  tfield2   = SourceQ_t_1d(x,t,A_x,A_t,D_t,k_0,cp_0,rho);

  tfield=tfield-tfield2;

  if(tfield > threshold)
    {
      cout << "\nMASA REGRESSION TEST FAILED: Heat Equation 1d Unsteady Constant\n";
      cout << "T Source Term\n";
      cout << "Exceeded Threshold by: " << tfield << endl;
      exit(1);
    }

  cout << "1D Unsteady Constant Heat Equation: PASSED\n";

  /// -----------------------------------------------------------------------
  // initalize 2D
  /// -----------------------------------------------------------------------
  masa_init("temp-test-2d","heateq_2d_unsteady_const");

  masa_set_param("A_x",param);
  A_x=param;

  masa_set_param("A_t",param);
  A_t=param;

  masa_set_param("D_t",param);
  D_t=param;

  masa_set_param("cp_0",param);
  cp_0=param;

  masa_set_param("rho",param);
  rho=param;

  masa_set_param("k_0",param);
  k_0=param;

  masa_set_param("B_y",param);
  double B_y=param;

  masa_set_param("B_t",param);
  double B_t=param;

  // evaluate source terms (2D)
  masa_sanity_check();
  masa_eval_t_source(x,y,t,&tfield);
  tfield2   = SourceQ_t_2d(x,y,t,A_x,A_t,B_y,B_t,D_t,rho,k_0,cp_0);

  tfield=tfield-tfield2;

  if(tfield > threshold)
    {
      cout << "\nMASA REGRESSION TEST FAILED: Heat Equation 2d Unsteady Constant\n";
      cout << "T Source Term\n";
      cout << "Exceeded Threshold by: " << tfield << endl;
      exit(1);
    }

  cout << "2D Unsteady Constant Heat Equation: PASSED\n";


  /// -----------------------------------------------------------------------
  // initalize 3D
  /// -----------------------------------------------------------------------
  masa_init("temp-test-3d","heateq_3d_unsteady_const");

  masa_set_param("A_x",param);
  A_x=param;

  masa_set_param("A_t",param);
  A_t=param;

  masa_set_param("D_t",param);
  D_t=param;

  masa_set_param("cp_0",param);
  cp_0=param;

  masa_set_param("rho",param);
  rho=param;

  masa_set_param("k_0",param);
  k_0=param;

  masa_set_param("B_y",param);
  B_y=param;

  masa_set_param("B_t",param);
  B_t=param;

  masa_set_param("C_z",param);
  double C_z=param;

  masa_set_param("C_t",param);
  double C_t=param;

  // evaluate source terms (3D)
  masa_sanity_check();
  masa_eval_t_source(x,y,z,t,&tfield);
  tfield2   = SourceQ_t_3d(x,y,z,t,A_x,A_t,B_y,B_t,C_z,C_t,D_t,k_0,cp_0,rho);

  tfield=tfield-tfield2;

  if(tfield > threshold)
    {
      cout << "\nMASA REGRESSION TEST FAILED: Heat Equation 3d Unsteady Constant\n";
      cout << "T Source Term\n";
      cout << "Exceeded Threshold by: " << tfield << endl;
      exit(1);
    }

  cout << "3D Unsteady Constant Heat Equation: PASSED\n";

  // presumably, all tests passed
  return 0;
}
