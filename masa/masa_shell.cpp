/*--------------------------------------------------------------------------
  *--------------------------------------------------------------------------
  *
  * Copyright (C) 2010 The PECOS Development Team
  *
  * Please see http://pecos.ices.utexas.edu for more information.
  *
  * This file is part of MASA.
  *
  * MASA is free software: you can redistribute it and/or modify it under
  * the terms of the GNU Lesser General Public License as published by the Free
  * Software Foundation, either version 3 of the License, or (at your option)
  * any later version.
  *
  * MASA is distributed in the hope that it will be useful, but WITHOUT ANY
  * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
  * details.
  *
  * You should have received a copy of the GNU Lesser General Public License along 
  * with MASA.  If not, see <http://www.gnu.org/licenses/>.
  *
  *--------------------------------------------------------------------------
  
  MASA -- Manufactured Analytical Solutions Abstraction Library

  A software interface that provides access to all manufactured solutions to 
  be used by various models throughout the center.
  
  *--------------------------------------------------------------------------
  */  

// this is the MASA shell, for users to grow familiar with the masa namespace, as well as debugging

#include "masa.h"

using namespace MASA;

void masa_shell_choose_solution()
{
  char userstring[100],masterstring[100];
  void* ptr; // initialize null
  double dbl=12.3;
  double dbl2;
  string str;

  printf("\n\n Now type solution to initialize.");
  printf("\n Remember to place underscores between parameters (e.g. steady_heatequation)");
  printf("\n\n");
  scanf("%s",userstring);

  //masa_map_dimension(userstring,masterstring);
  masa_getid(&ptr,userstring);
  masa_get_name(ptr,&str);
  cout << endl << "User has selected: " << str << endl;

  cout << "input double" << endl;
  cin >> dbl;
  cout << "double is:" << dbl << endl;
  masa_get_param(ptr,"axp",&dbl2);
  cout << "initially set to:" << dbl2 << endl;
  masa_set_param(ptr,"axp",dbl);
  masa_get_param(ptr,"axp",&dbl2);
  cout << "double is:" << dbl2 << endl;

}// end masa_shell_choose_solution

void masa_shell_print_avail()
{
  // this function prints the available manufactured solutions -- 
  // obviously needs to get a handle on the manufactured solutions class here
  masa_printid();
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
      scanf("%i",&q);
      
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


