// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// $Author$
// $Id$
//
// c_euler_example.c: this is an example of the heat equation in 1d for C
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <cmasa.h>
#include <stdio.h>

int main()
{
  double sol;
  
  // init
  cmasa_init("nick","euler_1d");
  cmasa_init_param();

  // list
  cmasa_list_mms();
  cmasa_display_param();

  //check all initialized properly
  cmasa_sanity_check();
  cmasa_eval_1d_u_source(1.2,&sol);
  printf("\nt source: %g\n",sol);

  return 0; // done
}
