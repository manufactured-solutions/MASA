#include <math.h>
#include <masa.h>

using namespace std;
using namespace MASA;

const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

double SourceQ_t_1d(
  double x,
  double A_x,
  double k_0,
  double k_1,
  double k_2)
{
  double Q_T = 0.3e1 * A_x * A_x * k_2 * pow(cos(A_x * x), 0.3e1) + 0.2e1 * A_x * A_x * k_1 * pow(cos(A_x * x), 0.2e1) - A_x * A_x * k_1 + (k_0 - 0.2e1 * k_2) * A_x * A_x * cos(A_x * x);
  return Q_T;
}

double SourceQ_t_2d (
  double x,
  double y,
  double A_x,
  double B_y,
  double k_0,
  double k_1,
  double k_2)
{
  double Q_T = (0.3e1 * A_x * A_x + 0.3e1 * B_y * B_y) * k_2 * pow(cos(A_x * x), 0.3e1) * pow(cos(B_y * y), 0.3e1) + (0.2e1 * A_x * A_x + 0.2e1 * B_y * B_y) * k_1 * pow(cos(A_x * x), 0.2e1) * pow(cos(B_y * y), 0.2e1) - (pow(cos(B_y * y), 0.2e1) * A_x * A_x + pow(cos(A_x * x), 0.2e1) * B_y * B_y) * k_1 + (k_0 * A_x * A_x + k_0 * B_y * B_y - 0.2e1 * pow(cos(B_y * y), 0.2e1) * k_2 * A_x * A_x - 0.2e1 * pow(cos(A_x * x), 0.2e1) * k_2 * B_y * B_y) * cos(A_x * x) * cos(B_y * y);
  return Q_T;
}

double SourceQ_t_3d (
  double x,
  double y,
  double z,
  double A_x,
  double B_y,
  double C_z,
  double k_0,
  double k_1,
  double k_2)
{
  double Q_T = (0.3e1 * A_x * A_x + 0.3e1 * B_y * B_y + 0.3e1 * C_z * C_z) * k_2 * pow(cos(A_x * x), 0.3e1) * pow(cos(B_y * y), 0.3e1) * pow(cos(C_z * z), 0.3e1) + (0.2e1 * A_x * A_x + 0.2e1 * B_y * B_y + 0.2e1 * C_z * C_z) * k_1 * pow(cos(A_x * x), 0.2e1) * pow(cos(B_y * y), 0.2e1) * pow(cos(C_z * z), 0.2e1) - (pow(cos(B_y * y), 0.2e1) * pow(cos(C_z * z), 0.2e1) * A_x * A_x + pow(cos(A_x * x), 0.2e1) * pow(cos(C_z * z), 0.2e1) * B_y * B_y + pow(cos(A_x * x), 0.2e1) * pow(cos(B_y * y), 0.2e1) * C_z * C_z) * k_1 + (k_0 * A_x * A_x + k_0 * B_y * B_y + k_0 * C_z * C_z - 0.2e1 * pow(cos(B_y * y), 0.2e1) * pow(cos(C_z * z), 0.2e1) * k_2 * A_x * A_x - 0.2e1 * pow(cos(A_x * x), 0.2e1) * pow(cos(C_z * z), 0.2e1) * k_2 * B_y * B_y - 0.2e1 * pow(cos(A_x * x), 0.2e1) * pow(cos(B_y * y), 0.2e1) * k_2 * C_z * C_z) * cos(A_x * x) * cos(B_y * y) * cos(C_z * z);
  return Q_T;
}

int main()
{
  double tfield;
  double tfield2;
  double param=1.2;
  double x=.5;
  double y=.4;
  double z=.3;

  /// -----------------------------------------------------------------------
  // initalize 1D
  /// -----------------------------------------------------------------------
  masa_init("temp-test-1d","heateq_1d_steady_var");

  masa_set_param("A_x",param);
  double A_x=param;

  masa_set_param("k_0",param);
  double k_0=param;

  masa_set_param("k_1",param);
  double k_1=param;

  masa_set_param("k_2",param);
  double k_2=param;

  // evaluate source terms (1D)
  masa_sanity_check();
  masa_eval_t_source(x,&tfield);
  tfield2   = SourceQ_t_1d(x,A_x,k_0,k_1,k_2);

  tfield=tfield-tfield2;

  if(tfield > threshold)
    {
      cout << "\nMASA REGRESSION TEST FAILED: Heat Equation 1d Steady Variable\n";
      cout << "ERROR in: T Source Term\n";
      cout << "Exceeded Threshold by: " << tfield << endl;
      exit(1);
    }

  cout << "1D Steady Variable Coefficient Heat Equation: PASSED\n";
  //cout << "Residual: "<< tfield << endl;

  /// -----------------------------------------------------------------------
  // initalize 2D
  /// -----------------------------------------------------------------------
  masa_init("temp-test-2d","heateq_2d_steady_var");

  masa_set_param("A_x",param);
  A_x=param;

  masa_set_param("B_y",param);
  double B_y=param;

  masa_set_param("k_0",param);
  k_0=param;

  masa_set_param("k_1",param);
  k_1=param;

  masa_set_param("k_2",param);
  k_2=param;

  // evaluate source terms (2D)
  masa_sanity_check();
  masa_eval_t_source(x,y,&tfield);
  tfield2   = SourceQ_t_2d(x,y,A_x,B_y,k_0,k_1,k_2);

  tfield=tfield-tfield2;

  if(tfield > threshold)
    {
      cout << "\nMASA REGRESSION TEST FAILED: Heat Equation 2d Steady Variable\n";
      cout << "ERROR in: T Source Term\n";
      cout << "Exceeded Threshold by: " << tfield << endl;
      exit(1);
    }

  cout << "2D Steady Variable Coefficient Heat Equation: PASSED\n";
  //cout << "Residual: "<< tfield << endl;

  /// -----------------------------------------------------------------------
  // initalize 3D
  /// -----------------------------------------------------------------------
  masa_init("temp-test-3d","heateq_3d_steady_var");

  masa_set_param("A_x",param);
  A_x=param;

  masa_set_param("B_y",param);
  B_y=param;

  masa_set_param("C_z",param);
  double C_z=param;

  masa_set_param("k_0",param);
  k_0=param;

  masa_set_param("k_1",param);
  k_1=param;

  masa_set_param("k_2",param);
  k_2=param;

  // evaluate source terms (2D)
  masa_sanity_check();
  masa_eval_t_source(x,y,z,&tfield);
  tfield2   = SourceQ_t_3d(x,y,z,A_x,B_y,C_z,k_0,k_1,k_2);

  tfield=tfield-tfield2;

  if(tfield > threshold)
    {
      cout << "\nMASA REGRESSION TEST FAILED: Heat Equation 3d Steady Variable\n";
      cout << "ERROR in: T Source Term\n";
      cout << "Exceeded Threshold by: " << tfield << endl;
      exit(1);
    }

  cout << "3D Steady Variable Coefficient Heat Equation: PASSED\n";
  //cout << "Residual: "<< tfield << endl;
}
