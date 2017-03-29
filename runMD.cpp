////////////////////////////////////////////////////
// 		CSE392 - project
// Primary control program for our MD simulator
//
// Authors: James Thorpe & Rohit Satija
// The Univserity of Texas at Austin, 2017
//
//
// Contact us at james.thorpe@utexas.edu
////////////////////////////////////////////////////

// predefined headers
#include <iostream>
#include <new>
#include <vector>
// our defined headers
#include "parser.hpp"
#include "killer.hpp"

int main()
{
  // Variable declarations
  // status 		: int, stores error and messages for output
  // N			: int, number of molecules in simulation (int)
  // sl			: side length of the box (float, nm) 
  // T			: temperature (float, K)
  // ts			: time step (float, ns)
  // ns			: number of time stimes (int)
  // options 		: 1D int array, stores extra options the user inputs
  
  // Variables for the simulation
  int status=0,N=10,ns=100;
  int options [1]={0};		
  float sl=10.0,T=298.15,ts=1.0;

  // Internal variables
  int i,j,k;

  // ~~~~~~~~~~~		Begin Program		~~~~~~~~~~//
  // Comments :

  std::cout <<  "Starting runMD, Version 0.0 ...." << std::endl;

  // Create our running objects
  Killer killer;

  // ~~~~~~~~~~			Get Input		~~~~~~~~~~//
  // Comments:

  Parser parser;	//this creates our "Parser" class object, "parser"
  status = parser.getInput(&N, &sl, &T, &ts, &ns, options);
  if (status != 0)
  {
    return status;
  }

  //Last line
  std::cout << std::endl << "Exiting runMD with status :" << status << std::endl;
  return status;

}

