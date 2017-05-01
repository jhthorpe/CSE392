// This initializes our simulation space
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <math.h>
#include <random>
#include "init.hpp"
using namespace std;

// tester function
void tester()
{
  cout << "foo?" << endl;

}

// Construct and populate the simulation space
// Currently, this creates an array of atoms
// Obviously, this will need to be improved later

int Init::initialize(int *N, double *sl, double *T,vector<double> *m, vector<double> *pos, vector<double> *vel)
{
  
  // Varaibles
  // N		: input number of molecules
  // sl		: input side length of box
  // T		: input temperature
  // m		: input mass
  // pos	: input 1D vector of positions
  // vel	: input 1D vector of velocities
  // bl		: calculated distances between atoms (overshoot)
  // pl		: atoms per line (overshoot)

  // Internal variables
  int i,flag;
  double kB = 1.38064852e-23; 		// in J/K, which coincidentely works for nm^2/ns^2

  cout << "Initializer called..." << endl;

  int pl = ceil(pow(*N,1.0/3.0));
  double bl = double(*sl /double(pl+1)); 

  cout << "particles will be placed " << bl << " nm apart." << endl;

  random_device generator;		//random_device class object, generator

  // One might argue that this is not the best approach for the velocities, but oh well
 // double stdv = pow(kB * *T / (*m)[0],0.5);
 // normal_distribution<double> vdist(0,stdv);
 
  // Create our cube of molecules 
  for (i=0; i < *N; i++) 
  {
    double stdv = pow(kB * *T / (*m)[i],0.5);
    normal_distribution<double> vdist(0,stdv);
    // I could, of course, have used 3 loops, but this is better flops/mops ;)

    (*pos)[3*i] = bl * double(i % pl);						// x position in nm 
    (*pos)[3*i+1] = bl * double((i / pl) % pl);					// y position in nm
    (*pos)[3*i+2] = bl * double((i / int(pow(pl,2.0))) % int(pow(pl,2.0)));	// z position in nm

    (*vel)[3*i]=vdist(generator);					// x velocity in nm/ns
    (*vel)[3*i+1]=vdist(generator);				// y velocity in nm/ns
    (*vel)[3*i+2]=vdist(generator);				// z velocity in nm/ns

  }

  // Write out initial position and velocities
  ofstream posfile; 		//ofstream class object, posfile, for initial positions
  ofstream velfile; 		//ofstream class object, velfile, for initial velocities

  posfile.open("init_pos.txt");
  velfile.open("init_vel.txt");

  posfile << "xpos, ypos, zpos (nm)\n";
  velfile << "x velocity, y velocity, z velocity (nm/ns)\n"; 
  for (i=0; i < *N ; i++)
  {
    posfile << (*pos)[3*i] << "  " << (*pos)[3*i+1] << "  " << (*pos)[3*i+2] << "\n";
    velfile << (*vel)[3*i] << "  " << (*vel)[3*i+1] << "  " << (*vel)[3*i+2] << "\n";
  }
  
  posfile.close();
  velfile.close();

  cout << "Initial positions written to init_pos.txt\nInitial velocities written to init_vel.txt" << endl;

  return 0;
}


