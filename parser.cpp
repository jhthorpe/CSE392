// Program to parser input from MD.dat
#include <iostream>
#include <fstream>
#include "parser.hpp"

int Parser::getInput(int *numMol, float boxdim[], int options[])
{
  std::cout << "Parser called. Getting input..." << std::endl;
  std::cout << "boxdim[0] = " << boxdim[0] << std::endl; 
  std::cout << "options[0] = " << options[0] << std::endl; 
  
  boxdim[1] = 42.05;
  options[0] = 22; 

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



