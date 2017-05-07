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
#include <mpi.h>
// our defined headers
#include "parser.hpp"
#include "killer.hpp"
#include "init.hpp"
#include "forces.hpp"
#include "forces_pvfmm.hpp"
#include "verlet.hpp"
using namespace std;


int main(int argc, char* argv[])
{
  // Variable declarations
  // status 		: int, stores error and messages for output
  // N			: int, number of molecules in simulation (int)
  // sl			: side length of the box (double, nm) 
  // T			: temperature (double, K)
  // ts			: time step (double, ns)
  // ns			: number of time stimes (int)
  // ss                 : number of time steps to display or write (int)
  // m			: double of input mass (double, kg)
  // kin                : total kinetic energy (double, kg-nm^2/ns^2)
  // LJpot              : total lennard jones potential energy (double, kg-nm^2/ns^2)
  // elcpot             : total electrostatic potential energy (double, kg-nm^2/ns^2)
  // mass		: 1D vector of mass (vector<doubles>) (kg/particle)  
  // pos		: 1D vector of positions, stored x(n),y(n+1),z(n+2) (vector<double>)
  // vel		: 1D vector of velocities, stored dx(n),dy(n+1),dz(n+2) (vector<double>, len/time)
  // q			: 1D vector of charges
  // options 		: 1D int array, stores extra options the user inputs
  
  // Variables for the simulation
  int status=0,N=10,ns=100,ss=10;
  int options [1]={0};		
  double sl=10.0,T=298.15,ts=1.0,sig=1.0,eps=0.0,m=1.0,chrg=0.0,kin=0.0,LJpot=0.0,elcpot=0.0;

  // Internal variables
  int i,j,k;

  // ~~~~~~~~~~~		Begin Program		~~~~~~~~~~//
  // Comments :

  cout <<  "Starting runMD, Version 0.0 ...." << endl;

  // Create our running objects
  Killer killer;
  Init init;

  // ~~~~~~~~~~			Get Input		~~~~~~~~~~//
  // Comments: May want to test that the types are correct - Mar 28, 2017

  Parser parser;	//this creates our "Parser" class object, "parser"
  status = parser.getInput(&N, &sl, &T, &m, &ts, &ns, &ss, &sig, &eps, &chrg, options);
  if (status != 0)
  {
    killer.kill(status);
    return status;
  }

  // ~~~~~~~~~~			Initialize Box		~~~~~~~~~~//
  // Comments: Needs to have parallel treatment. Velocities are in nm^2/ns^2. 
  // I have essentially hardcoded
  i = N * 3;
  vector<double> pos;
  vector<double> vel;
  vector<double> mass(N, m);		//vector of masses, not efficient or flexible right now
  pos.reserve(N*3);	//reserve, but do not initialize, N*3 space. More efficient.
  vel.reserve(N*3);
  

  Init builder;

  status = builder.initialize(&N,&sl,&T,&mass,&pos,&vel);
  if (status != 0)
  {
    killer.kill(status);
  }

  // ~~~~~~~~~~			Testing forces		~~~~~~~~~~//
  // Comments: make sure the directionality is being handled correctly... also add in boundary conditions 
  // This currently is NOT set up for boudnary conditions, the first particle is set at 0,0. Fix this.

  MPI_Init(&argc, &argv);
  MPI_Comm comm=MPI_COMM_WORLD;

  vector<double> force(N*3, 0.0);
  
  Forces_pvfmm forces_pvfmm;	//Forces class object, forces

  //forces.LJ_seq_bound(&N,&sl,&sig,&eps,&pos,&force);
  //cout << "seq done\n";
//  forces.LJ_omp_bound(&N,&sl,&sig,&eps,&pos,&force);
  
  vector<double> q(N, chrg);
  forces_pvfmm.elc_pvfmm(&N,&sl,&q,&pos,&force,&comm);

  // ~~~~~~~~~~			Start verlet		~~~~~~~~~~//
  // Comments: sequential velocity verlet integration. integrates position and velocity with periodic boundary conditions

//  verlet v;           //verlet class object, v

//  v.Integration(&N,&ns,&ss,&ts,&sl,&sig,&eps,&q,&mass,&pos,&vel);

  //energy e;          //energy class object, e

  //e.kinetic_seq(&N,&mass,&vel,&kin,&T);
  //e.LJpot_seq(&N,&sl,&sig,&eps,&pos,&LJpot);
  //e.elcpot_seq(&N,&q,&pos,&elcpot);

  //Last line
  cout << endl << "Exiting runMD with status :" << status << endl;

  MPI_Finalize();

  return status;

}

