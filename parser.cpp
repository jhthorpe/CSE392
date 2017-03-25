// Program to parser input from MD.dat
#include <iostream>
#include <fstream>
#include "parser.hpp"
using namespace std;

int Parser::getInput(int *numMol, float boxdim[], int options[])
{
  std::cout << "Parser called. Getting input..." << std::endl;

  // Open input.dat
  fstream inFile;	// Create fstream class object, "inFile"
  inFile.open ("input.dat", ios::in);
  
  // Check the file was opened correctly
  if (! inFile.is_open() )
  {
    cout << "Parser could not open input.dat. Exiting." << endl;
    return 1;
  }

  

  // Close input.dat  
  inFile.close();

  return 0; 
} 

// my tester funtion, will later do stuff 
int Parser::parse(int *a, int *b)
{
  std::cout << "hi there" << std::endl;
  *a = 5;
  *b = 20;
  return 0; 
}



