// Program to parser input from MD.dat
#include <iostream>
#include <fstream>
#include "parser.hpp"
using namespace std;

int Parser::getInput(int *N, double *sl, double *T, double *m, double *ts, int *ns, double *sig, double *eps, int options[])
{
  int i;

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
    inFile >> *N >> *sl >> *T >> *m >> *ts >> *ns >> *sig >> *eps; 
  }

  // Close input.dat  
  inFile.close();

  // Print what saw
  cout << "Parser saw these values:" << endl;
  cout << "Number of molecules : " << *N << endl;
  cout << "Box length (nm) : " << *sl << endl; 
  cout << "Temperature (K) : " << *T << endl; 
  cout << "Molecule mass (g/mol) : " << *m << endl; 
  cout << "Time Step (fs) : " << *ts << endl; 
  cout << "Number of Steps : " << *ns << endl; 
  cout << "LJ sigma (m) : " << *sig << endl;
  cout << "LJ epsilon (J) : " << *eps << endl;

  // Correct units
  cout << endl;
  cout << "Converting units..." << endl;
  cout <<  "Mass -> kg/particle" << endl;
  *m = *m / 6.02214e26; 
  cout << "new mass: " << *m << endl;

  cout << "Meters -> nm" << endl;
  *sig = *sig * 1.0e9;
  cout << "new LJ sigma : " << *sig << endl; 


  cout << "====================" << endl; 
   
  return 0; 
} 


