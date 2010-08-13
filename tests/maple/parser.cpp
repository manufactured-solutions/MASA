#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <stdlib.h>
#include <sstream>

using namespace std;

int main(int argc, char* argv[])
{
  int p,p2,dimension,numparam;
  char ch;
  string str,substr,substr2,substr3,name;
  string comma = ","; 
  //stringstream ss;

  // need one input argument
  if(argc != 2) 
    {
      cout << "wrong number of arguments\n";
      exit(1);
    }

  ifstream from(argv[1]);
  if(!from) 
    {
      cout << "cannot open file!\n";
      exit(1);
    }

  // get first line
  getline(from,str);
  p = str.find(comma); //find first comma
  stringstream ss(str.substr(0,p));
  ss >> dimension;    
  cout << dimension << endl;

  p2 = str.find(comma,p+1); // 2nd comma
  name = str.substr(p+1,p2-p-1);
  cout << name << endl;
  string start = str.substr(0, 5);  

  cout << str.substr(p2+1) << endl;
  ss << str.substr(p2+1,1); // final chunk
  ss >> numparam;
  cout << numparam << endl;
  
  if(!from.eof()) cout << "something very weird has happened: file did not end!";
  
}
