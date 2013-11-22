// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012,2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
// $Author$
// $Id$
//
// masa_shell.cpp: for users to grow familiar with the masa namespace
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa.h>
#include <cstdio> 
#include <iostream>
#include <string>

using namespace std;
using namespace MASA;

typedef double Scalar;

void MASA::masa_shell_choose_solution()
{
  char userstring[100];

  int q = 1;
  int dimension;

  Scalar x,y,z,t;
  Scalar dbl;
  Scalar dbl2;
  Scalar field;

  string str1;
  
  printf("\n\n Now type solution to initialize.");
  printf("\n Remember to place underscores between parameters (e.g. heateq_1d_steady_const)");
  printf("\n\n");
  cin >> userstring;

  //masa_map_dimension(userstring,masterstring);
  masa_init<Scalar>("dummy",userstring);
  masa_get_name<Scalar>(&str1);
  cout << endl << "User has selected: " << str1 << endl;

  // now let user register variables, etc.

  while(q>0)
    {
      printf("\n\n 0) Exit and drop currently selected manufactured solution");
      printf("\n 1) Show all variables and registered values");
      printf("\n 2) Register variable values");
      printf("\n 3) Set to default values");
      printf("\n 4) Evaluate\n\n");
      cin >> q;

      // switch statement based on user input
      switch(q)
	{
	case 0:
	  printf("\n Exiting.");
	  break;

	case 1:
	  printf("\n User Selected 1: Display all variables\n");
	  masa_display_param<Scalar>();
	  break;
	  
	case 2:
	  printf("\n User Selected 2: Register Variable");
	  printf("\n Input variable name:\n");
	  cin >> userstring;
	  dbl2 = masa_get_param<Scalar>(userstring);
	  cout << "currently set to:" << dbl2 << endl;
	  cout << "\nInput new value (double)" << endl;
	  cin >> dbl;
	  masa_set_param<Scalar>(userstring,dbl);
	  dbl2 = masa_get_param<Scalar>(userstring);
	  cout << endl << userstring << " is now set to:" << dbl2 << endl;
	  break;

	case 3:
	  printf("\n User Selected 3: Set to default values\n");
	  break;
	  
	case 4:
	  printf("\n User Selected 4: Evaluate");
	  masa_sanity_check<Scalar>();
	  masa_get_dimension<Scalar>(&dimension);	 
	  switch(dimension)
	    {
	    case 1:
	      cout << "\nplease input x location: \n";
	      cin >> x;
	      
	      field = masa_eval_source_t<Scalar>(x);
	      cout << "source term is:" << field;
	      break;

	    case 2:
	      cout << "\nplease input x location: \n";
	      cin >> x;

	      cout << "\nplease input y location: \n";
	      cin >> y;

	      field = masa_eval_source_u<Scalar>(x,y);
	      cout << "source term is:" << field;
	      break;

	    case 3:
	      cout << "\nplease input x location: \n";
	      cin >> x;

	      cout << "\nplease input y location: \n";
	      cin >> y;

	      cout << "\nplease input z location: \n";
	      cin >> z;

	      field = masa_eval_source_u<Scalar>(x,y,z);
	      cout << "source term is:" << field << endl << endl;
	      break;
	      
	    case 4:
	      cout << "\nplease input x location: \n";
	      cin >> x;

	      cout << "\nplease input y location: \n";
	      cin >> y;

	      cout << "\nplease input z location: \n";
	      cin >> z;

	      cout << "\nplease input t location: \n";
	      cin >> t;

	      cout << "source term rho   is:" << masa_eval_source_rho  <Scalar>(x,y,z,t) << endl;
	      cout << "source term rho_u is:" << masa_eval_source_rho_u<Scalar>(x,y,z,t) << endl;
	      cout << "source term rho_v is:" << masa_eval_source_rho_v<Scalar>(x,y,z,t) << endl;
	      cout << "source term rho_w is:" << masa_eval_source_rho_w<Scalar>(x,y,z,t) << endl;
	      cout << "source term rho_e is:" << masa_eval_source_rho_e<Scalar>(x,y,z,t) << endl;
	      cout << endl;
	      break;
	      
	    default: 
	      printf("\n Error: Undefined dimension for class, please try again.\n");
	      break;

	    }//done with switch
	  break; // done with case 3
	  
	default:
	  printf("\n Error: Undefined input, please try again.\n");
	  break;
	} // end switch
     
    }

}// end masa_shell_choose_solution

void MASA::masa_shell_print_avail()
{
  // this function prints the available manufactured solutions -- 
  // obviously needs to get a handle on the manufactured solutions class here
  masa_printid<Scalar>();
}

int main()
{
  int q; // this defines the while loop

  // welcome message
  printf("\n\n* ----------------------------------------------------------------------------- *");
  printf("\n* Welcome to the MASA (Manufactured Analytical Solutions Abstraction) Library.");
  printf("\n* This is a software shell interface which will provide access to all");
  printf("\n* manufactured solutions provided in this library.");
  printf("\n* ----------------------------------------------------------------------------- *\n\n");

  q=1600;

  //enter shell
  while(q > 0)
    {
      
      // user chooses option

      printf("\n 1) Display Available Solutions");
      printf("\n 2) Choose and Initialize Solution");
      printf("\n Select option (0 to exit)");
      printf("\n input: ");
      cin >> q;
      
      // switch statement based on user input
      switch(q)
	{
	case 0:
	  printf("\n Exiting.");
	  break;

	case 1:
	  printf("\n User Selected 1: Display Available Solutions");
	  masa_shell_print_avail();
	  break;
	  
	case 2:
	  printf("\n User Selected 2: Choose and Initialize Solution");
	  masa_shell_choose_solution();
	  break;
	  
	default:
	  printf("\n Error: Undefined input, please try again.\n");
	  break;
	} // end switch
      
    } // end shell

  printf("\n Thank you for using the MASA shell. Have a nice day.");
  printf("\n* ----------------------------------------------------------------------------- *\n\n");
} // end main


