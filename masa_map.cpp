// this routine takes a string, parses it
// and returns what we 'really' wanted the user to input

#include "masa.h"

// user chooses 1d,2d, or 3d 
int masa_map_dimension(string input_string, string return_string)
{
  multimap<string, string> names;
  int ret_val=0;  // exits with zero upon failure, 1 for success
  //map<string, string> names;

  // map 1-d
  names.insert(pair<string, string>("1D",               "1d"));
  names.insert(pair<string, string>("1d",               "1d"));
  names.insert(pair<string, string>("one dimensional",  "1d"));
  names.insert(pair<string, string>("one-dimensional",  "1d"));
  names.insert(pair<string, string>("one",              "1d"));
  names.insert(pair<string, string>("oned",             "1d"));
  names.insert(pair<string, string>("one-d",            "1d"));
  names.insert(pair<string, string>("one-D",            "1d"));
  names.insert(pair<string, string>("oneD",             "1d"));

  // map 2-d
  names.insert(pair<string, string>("2D",               "2d"));
  names.insert(pair<string, string>("2d",               "2d"));
  names.insert(pair<string, string>("two dimensional",  "2d"));
  names.insert(pair<string, string>("two-dimensional",  "2d"));
  names.insert(pair<string, string>("two-dim",          "2d"));
  names.insert(pair<string, string>("two",              "2d"));
  names.insert(pair<string, string>("twod",             "2d"));
  names.insert(pair<string, string>("two-d",            "2d"));
  names.insert(pair<string, string>("two-D",            "2d"));
  names.insert(pair<string, string>("twoD",             "2d"));

  // map 3-d
  names.insert(pair<string, string>("3D",               "3d"));
  names.insert(pair<string, string>("3d",               "3d"));
  names.insert(pair<string, string>("three dimensional","3d"));
  names.insert(pair<string, string>("three-dimensional","3d"));
  names.insert(pair<string, string>("three-dim",        "3d"));
  names.insert(pair<string, string>("three",            "3d"));
  names.insert(pair<string, string>("threed",           "3d"));
  names.insert(pair<string, string>("three-d",          "3d"));
  names.insert(pair<string, string>("three-D",          "3d"));
  names.insert(pair<string, string>("threeD",           "3d"));

  multimap<string, string>::iterator p;

  p = names.find(input_string);
  if(p != names.end()) { // found a name
    do {
      cout << input_string << ", " << p->second;
      cout << endl;
      p++;
    } while (p != names.upper_bound(input_string));
    ret_val=1;
  }
  else{
    cout << "Name not found.\n";
    ret_val = 0;
  }

  return ret_val;

}// end masa_map_dimension
