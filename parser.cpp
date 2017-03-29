// Program to parser input from MD.dat
#include <iostream>
#include <fstream>
#include "parser.hpp"
using namespace std;

int Parser::getInput(int *N, float *sl, float *T, float *ts, int *ns, int options[])
{
  cout << "====================" << endl; 
  cout << "Parser called. Getting input..." << endl;

  // Open input.dat
  fstream inFile;	// Create fstream class object, "inFile"
  inFile.open ("input.dat", ios::in);
  
  // Check the file was opened correctly
  if (! inFile.is_open() )
  {
    cout << "Parser could not open input.dat. Exiting." << endl;
    return 1;
  }
 
  // otherwise, read in values 
  else
  {
    // Get initial values
    inFile >> *N >> *sl >> *T >> *ts >> *ns; 
  }

  // Close input.dat  
  inFile.close();

  // Print what saw
  cout << "Parser saw these values:" << endl;
  cout << "Number of molecules : " << *N << endl;
  cout << "Box length (nm) : " << *sl << endl; 
  cout << "Temperature (K) : " << *T << endl; 
  cout << "Time Step (fs) : " << *ts << endl; 
  cout << "Number of Steps : " << *ns << endl; 

  cout << "====================" << endl; 
   
  return 0; 
} 

