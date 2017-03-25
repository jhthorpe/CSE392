// Program that handles killing the rest of runMD 
// 	This is not currently anything fancy, but I may want it to be later
#include <iostream>
using namespace std;

int Killer::kill (int status)
{
 cout << "runMD killed with code :" << status << endl;
 return 0;
}
