
// this is the MASA shell, for users to grow familiar with the masa namespace, as well as debugging
// this probably should not be contained in the MASA namespace

#include "masa.h"

using namespace std;

void masa_shell_choose_solution()
{
  char userstring[100],masterstring[100];

  printf("\n\n Now type solution to initialize.");
  printf("\n Remember to place underscores between parameters (e.g. steady_heatequation)");
  printf("\n\n");
  scanf("%s",userstring);
  masa_map_dimension(userstring,masterstring);  
}// end masa_shell_choose_solution

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


