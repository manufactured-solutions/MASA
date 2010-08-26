// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// $Author$
// $Id$
//
// switch.cpp:
// C++ example demonstrating how to switch between two different solutions
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

// header contains all subroutine information
#include <masa.h>

// MASA was designed with a custom namespace
using namespace MASA;
using namespace std;

int main()
{
  double solution;

  // print all solutions
  masa_printid();
  cout << endl << endl;

  // initialize first solution
  masa_init("alice","heateq_1d_steady_const");
  masa_init_param();


  // initialize 2nd solution
  masa_init("bob"  ,"euler_2d");
  masa_init_param();
  
  // list all initialized solutions
  masa_list_mms();

  // lets manipulate the alice set
  masa_select_mms("alice");
  masa_display_param();
  masa_eval_t_source(1.2,&solution);
  cout << solution << endl;

  // now switch to and edit bob
  masa_select_mms("bob");
  masa_display_param();
  masa_eval_u_source(1,1,&solution);
  cout << solution << endl;

}

