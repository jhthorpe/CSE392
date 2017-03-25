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
  // numMol		: int, number of molecules in simulation
  // boxdim		: 1D float array, dimensions of our box
  // options 		: 1D int array, stores extra options the user inputs
  
  // Variables for the simulation
  int status=0,numMol=0;
  int options [1]={0};		
  float boxdim [3]={0.0,0.0,0.0};

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
  status = parser.getInput(&numMol, boxdim, options);
  if (status != 0)
  {
    return status;
  }
  
  std::cout << "boxdim[1] = " << boxdim[1] << std::endl;
  std::cout << "options[0] = " << options[0] << std::endl;
  

  //Last line
  std::cout << std::endl << "Exiting runMD with status :" << status << std::endl;
  return status;

}

