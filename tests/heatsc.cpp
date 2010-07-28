#include <math.h>
#include <masa.h>

using namespace std;
using namespace MASA;

const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

double SourceQ_t_1d(double x, double A_x, double k_0)
{
  double Q_T = A_x * A_x * k_0 * cos(A_x * x);
  return Q_T;
}

double SourceQ_t_2d (
  double x,
  double y,
  double A_x,
  double B_y,
  double k_0)
{
  double Q_T = k_0 * cos(A_x * x) * cos(B_y * y) * (A_x * A_x + B_y * B_y);
  return Q_T;
}

double SourceQ_t_3d (
  double x,
  double y,
  double z,
  double A_x,
  double B_y,
  double C_z,
  double k_0)
{
  double Q_T = k_0 * cos(A_x * x) * cos(B_y * y) * cos(C_z * z) * (A_x * A_x + B_y * B_y + C_z * C_z);
  return Q_T;
}

int main()
{
  double tfield,tfield2;
  double param=1.2;
  double x=.5;
  double y=.4;
  double z=.3;

  /// -----------------------------------------------------------------------
  // initalize 1D
  /// -----------------------------------------------------------------------
  masa_init("temp-test-1d","heateq_1d_steady_const");

  masa_set_param("A_x",param);
  double A_x=param;

  masa_set_param("k_0",param);
  double k_0=param;

  // evaluate source terms (1D)
  masa_sanity_check();
  masa_eval_t_source(x,&tfield);
  tfield2   = SourceQ_t_1d(x,A_x,k_0);

  tfield=tfield-tfield2;

  if(tfield > threshold)
    {
      cout << "\nMASA REGRESSION TEST FAILED: Heat Equation 1d Steady Constant\n";
      cout << "T Source Term\n";
      cout << "Exceeded Threshold by: " << tfield << endl;
      exit(1);
    }

  cout << "1D Steady Constant Heat Equation: PASSED\n";


  /// -----------------------------------------------------------------------
  // initalize 2D
  /// -----------------------------------------------------------------------
  masa_init("temp-test-2d","heateq_2d_steady_const");

  masa_set_param("A_x",param);
  A_x=param; // A_x already init

  masa_set_param("B_y",param);
  double B_y=param;

  masa_set_param("k_0",param);
  k_0=param; // k_0 already init

  // evaluate source terms (1D)
  masa_sanity_check();
  masa_eval_t_source(x,y,&tfield);
  tfield2   = SourceQ_t_2d(x,y,A_x,B_y,k_0);

  tfield=tfield-tfield2;

  if(tfield > threshold)
    {
      cout << "\nMASA REGRESSION TEST FAILED: Heat Equation 2d Steady Constant\n";
      cout << "T Source Term\n";
      cout << "Exceeded Threshold by: " << tfield << endl;
      exit(1);
    }

  cout << "2D Steady Constant Heat Equation: PASSED\n";

  /// -----------------------------------------------------------------------
  // initalize 3D
  /// -----------------------------------------------------------------------

  masa_init("temp-test-3d","heateq_3d_steady_const");

  masa_set_param("A_x",param);
  A_x=param; // A_x already init

  masa_set_param("B_y",param);
  B_y=param; // already init

  masa_set_param("C_z",param);
  double C_z=param;

  masa_set_param("k_0",param);
  k_0=param; // k_0 already init

  // evaluate source terms (1D)
  masa_sanity_check();
  masa_eval_t_source(x,y,z,&tfield);
  tfield2   = SourceQ_t_3d(x,y,z,A_x,B_y,C_z,k_0);

  tfield=tfield-tfield2;

  if(tfield > threshold)
    {
      cout << "\nMASA REGRESSION TEST FAILED: Heat Equation 3d Steady Constant\n";
      cout << "T Source Term\n";
      cout << "Exceeded Threshold by: " << tfield << endl;
      exit(1);
    }

  cout << "3D Steady Constant Heat Equation: PASSED\n";

  // presumably, all tests passed
  return 0;

}
