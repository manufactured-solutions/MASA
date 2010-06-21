#include <cmasa.h>
#include <masa.h>
#include <string>

using namespace std;
using namespace MASA;

char *asd(char* in, char *out)
{
  strcat(out, in); // <-- err arg 2 makes pointer from integer without a cast
  return out;
}

int cmasa_init(char* specificname,char* functionname)
{
  string sn(specificname);
  string fn(functionname);

  masa_init(sn,fn);

  return 0;
}

int cmasa_select_mms(char* function_user_wants)
{
  string fuw(function_user_wants);
  masa_select_mms(fuw);
  return 0;
}

/*
int MASA::cmasa_curr_mms(char* function_name)
{
  string fn(function_name);
  masa_curr_mms(fn);
  return 0;
}
*/

int cmasa_list_mms()
{
  masa_list_mms();
  return 0;
}

int cmasa_init_param()
{
  masa_init_param();
  return 0;
}

int cmasa_sanity_check()
{
  masa_sanity_check();
  return 0;
}

int cmasa_display_param()
{
  masa_display_param();
  return 0;
}
