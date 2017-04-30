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
#include <cmath>
// our defined headers
#include "parser.hpp"
#include "killer.hpp"
#include "init.hpp"
#include "forces.hpp"
using namespace std;

int main()
{
  // Variable declarations
  // status 		: int, stores error and messages for output
  // N			: int, number of molecules in simulation (int)
  // sl			: side length of the box (float, nm) 
  // T			: temperature (float, K)
  // ts			: time step (float, ns)
  // ns			: number of time stimes (int)
  // m			: mass (kg/particle) 
  // pos		: 1D vector of positions, stored x(n),y(n+1),z(n+2) (vector<float>)
  // vel		: 1D vector of velocities, stored dx(n),dy(n+1),dz(n+2) (vector<float>, len/time)
  // options 		: 1D int array, stores extra options the user inputs
  
  // Variables for the simulation
  int status=0,N=10,ns=100;
  int options [1]={0};		
  float sl=10.0,m=1.0,T=298.15,ts=1.0,sig=1.0,eps=0.0;

  // Internal variables
  int i,j,k;

  // ~~~~~~~~~~~		Begin Program		~~~~~~~~~~//
  // Comments :

  std::cout <<  "Starting runMD, Version 0.0 ...." << std::endl;

  // Create our running objects
  Killer killer;
  Init init;

  // ~~~~~~~~~~			Get Input		~~~~~~~~~~//
  // Comments: May want to test that the types are correct - Mar 28, 2017

  Parser parser;	//this creates our "Parser" class object, "parser"
  status = parser.getInput(&N, &sl, &T, &m, &ts, &ns, &sig, &eps, options);
  if (status != 0)
  {
    killer.kill(status);
    return status;
  }

  // ~~~~~~~~~~			Initialize Box		~~~~~~~~~~//
  // Comments: Needs to have parallel treatment. Velocities are in nm^2/ns^2. 
  i = N * 3;
  vector<float> pos;
  vector<float> vel;
  pos.reserve(N*3);
  vel.reserve(N*3);

  Init builder;

  status = builder.initialize(&N,&sl,&T,&m,&pos,&vel);
  if (status != 0)
  {
    killer.kill(status);
  }

  // ~~~~~~~~~~			Testing forces		~~~~~~~~~~//
  // Comments: 

  vector<float> force;
  force.reserve(N*3);
  
  Forces forces;	//Forces class object, forces

  forces.Test();
  forces.LJ_seq(&N,&sl,&sig,&eps,&pos,&force);

  // ~~~~~~~~~~			Start verlet		~~~~~~~~~~//
  // Comments:

  //Last line
  std::cout << std::endl << "Exiting runMD with status :" << status << std::endl;
  return status;

}

