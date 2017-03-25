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

int main()
{
  // Variable declarations
  // status 		: int, stores error and messages for output
  // numMol		: int, number of molecules in simulation
  // boxdim		: 1D float array, dimensions of our box
  // options 		: 1D int array, stores extra options the user inputs
  
  // Variables for the simulation
  int status,numMol;
  int options [1]={0};		
  float boxdim [3]={0.0,0.0,0.0};

  // Internal variables
  int i,j,k;

  // ~~~~~~~~~~~		Begin Program		~~~~~~~~~~//
  //First line
  std::cout <<  "Starting runMD, Version 0.0 ...." << std::endl;


  Parser parser;	//this creates our "Parser" class object, "parser"
  
  parser.getInput(&numMol, boxdim[], options[]);
  


  // Testing zone ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int a = 0, b = 0;
  std::cout << a << ", " << b << std::endl;
  parser.parse(&a,&b);	// pass the POINTERS &a and &b so we can act on their address 
  std::cout << a << ", " << b << std::endl;
  // END zone ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //Last line
  std::cout << std::endl << "Exiting runMD with status :" << status << std::endl;

}

