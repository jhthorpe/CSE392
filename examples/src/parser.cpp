// Program to parser input from MD.dat
#include <iostream>
#include <fstream>
#include <cmath>
#include <mpi.h>
#include "parser.hpp"
using namespace std;

int Parser::getInput(int *N, double *sl, double *T, double *m, double *ts, int *ns, int *ss, double *sig, double *eps, double *q, MPI_Comm * comm)
{
  int i,rank;

  MPI_Comm_rank(*comm,&rank);
  if (!rank) {
    cout << "====================" << endl; 
    cout << "Parser called. Getting input..." << endl;
  }

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
    inFile >> *N >> *sl >> *T >> *m >> *ts >> *ns >> *ss >> *sig >> *eps >> *q; 
  }

  // Close input.dat  
  inFile.close();

  
  // Print what saw
  if (!rank) {
    cout << "Parser saw these values:" << endl;
    cout << "Number of molecules : " << *N << endl;
    cout << "Box length (nm) : " << *sl << endl; 
    cout << "Temperature (K) : " << *T << endl; 
    cout << "Molecule mass (g/mol) : " << *m << endl; 
    cout << "Time Step (fs) : " << *ts << endl; 
    cout << "Number of Steps : " << *ns << endl; 
    cout << "Numner of steps to display : " << *ss << endl;
    cout << "LJ sigma (m) : " << *sig << endl;
    cout << "LJ epsilon (J) : " << *eps << endl;
    cout << "Charges : " << *q << endl;

  // Correct units
  cout << endl;
  cout << "Converting units..." << endl;
  cout <<  "Mass -> kg/particle" << endl;
  }
  *m = *m / 6.02214e26; 
  if (!rank) {
    cout << "new mass: " << *m << endl;

    cout << "Meters -> nm" << endl;
  }
  *sig = *sig * 1.0e9;
  if (!rank) {
    cout << "new LJ sigma : " << *sig << endl; 

    cout << "Charges -> Coulombs " << endl;
  }
  *q = *q * 1.9e-19;
  if (!rank) {
    cout << "new charges : " << *q << endl;

    cout << "fs -> ns " << endl;
  }
  *ts = *ts * 1e-6;
  if (!rank) {
    cout << "new times : " << *ts << endl;

  //Normalizing to side length
    cout << "Normalizing to side length. Following values will be changed:\n";
  }
  *sig = *sig / *sl; 
  *eps = *eps / pow(*sl,2);
  if (!rank) {
    cout << "LJ sigma, LJ epsilon, ke, position, velocity, force, energy ..." << endl; 
  
  
    cout << "====================" << endl; 
  }
   
  return 0; 
} 


