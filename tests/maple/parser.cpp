#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>

using namespace std;

int main(int argc, char* argv[])
{
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

  char ch;
  string str;
  //while(from.get(ch)) cout << ch;
  while(!from.eof())
    {
      getline(from,str);
      cout << str;
	      //cout << ch;
    }
  
  if(!from.eof()) cout << "something very weird has happened";
  
}
