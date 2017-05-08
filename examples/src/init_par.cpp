// This initializes our simulation space
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <math.h>
#include <random>
#include <mpi.h>
#include "init_par.hpp"
using namespace std;

// Construct and populate the simulation space
// Currently, this creates an array of atoms
// Obviously, this will need to be improved later

int Init_par::initialize_mpi(int *N, double *sl, double *T,vector<double> *m, vector<double> *pos, vector<double> *vel, int *start, int *end, MPI_Comm *comm)
{
  
  // Varaibles
  // N		: input number of molecules in this task
  // sl		: input side length of box
  // T		: input temperature
  // m		: input mass
  // pos	: input 1D vector of positions
  // vel	: input 1D vector of velocities
  // bl		: calculated distances between atoms (overshoot)
  // pl		: atoms per line (overshoot)
  // start	: starting index for task
  // end	: end index for task

  // Internal variables
  int i,flag;
  double kB = 1.38064852e-23; 		// in J/K, which coincidentely works for nm^2/ns^2
  int rank;

  MPI_Comm_rank(*comm,&rank);

  cout << "Initializer called..." << endl;

  int pl = ceil(pow(*N,1.0/3.0));
  double bl = double((*sl-0.1) /double(pl+1)); 

  cout << "particles will be placed " << bl << " nm apart." << endl;

  random_device generator;		//random_device class object, generator

  // Create our cube of molecules 
  for (i=*start; i < *end; i++) 
  {
    double stdv = pow(kB * *T / (*m)[i],0.5);
    normal_distribution<double> vdist(0,stdv);
    // I could, of course, have used 3 loops, but this is better flops/mops ;)

    (*pos)[3*(i-*start)] = (0.1 + bl * double(i % pl)) / *sl;						// x position in nm/sl 
    (*pos)[3*(i-*start)+1] = (0.1 + bl * double((i / pl) % pl)) / *sl;					// y position in nm/sl
    (*pos)[3*(i-*start)+2] = (0.1 +bl * double((i / int(pow(pl,2.0))) % int(pow(pl,2.0)))) / *sl;	// z position in nm/sl

    (*vel)[3*(i-*start)]=(vdist(generator)) / *sl;					// x velocity in nm/ns
    (*vel)[3*(i-*start)+1]=(vdist(generator)) / *sl;				// y velocity in nm/ns
    (*vel)[3*(i-*start)+2]=(vdist(generator)) / *sl;				// z velocity in nm/ns

  }

  // Write out initial position and velocities

  if (rank == 0 || rank == 1){

    ofstream posfile; 		//ofstream class object, posfile, for initial positions
    ofstream velfile; 		//ofstream class object, velfile, for initial velocities

    if (!rank){
      posfile.open("init_pos0.txt");
      velfile.open("init_vel0.txt");
    } else {
      posfile.open("init_pos1.txt");
      velfile.open("init_vel1.txt");
    }

    posfile << "xpos, ypos, zpos (nm)\n";
    velfile << "x velocity, y velocity, z velocity (nm/ns)\n"; 
    for (i=0; i < *N ; i++)
    {
      posfile << (*pos)[3*i]* *sl << "  " << (*pos)[3*i+1]* *sl << "  " << (*pos)[3*i+2]* *sl << "\n";
      velfile << (*vel)[3*i]* *sl << "  " << (*vel)[3*i+1]* *sl << "  " << (*vel)[3*i+2]* *sl << "\n";
    }
  
    posfile.close();
    velfile.close();

    cout << "Initial positions written to init_pos.txt\nInitial velocities written to init_vel.txt" << endl;

  }

  return 0;
}


