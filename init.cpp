// This initializes our simulation space
#include <iostream>
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

int Init::initialize(int *N, float *sl, float *T,float *m, vector<float> *pos, vector<float> *vel)
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
  float kB = 1.38064852e-23; 		// in J/K, which coincidentely works for nm^2/ns^2

  cout << "Initializer called..." << endl;

  int pl = ceil(pow(*N,1.0/3.0));
  float bl = float(*sl /float(pl+1)); 

  // One might argue that this is not the best approach for the velocities, but oh well
  float stdv = pow(kB * *T / *m,0.5);

  random_device generator;
  normal_distribution<float> vdist(0,stdv);
 
  // Create our cube of molecules 
  for (i=0; i < *N; i++) 
  {
    // I could, of course, have used 3 loops, but this is better flops/mops ;)

    (*pos)[3*i] = bl * float(i % pl);						// x position in nm 
    (*pos)[3*i+1] = bl * float((i / pl) % pl);					// y position in nm
    (*pos)[3*i+2] = bl * float((i / int(pow(pl,2.0))) % int(pow(pl,2.0)));	// z position in nm

    (*vel)[3*i]=vdist(generator);					// x velocity in nm/ns
    (*vel)[3*i+1]=vdist(generator);				// y velocity in nm/ns
    (*vel)[3*i+2]=vdist(generator);				// z velocity in nm/ns

  }

  return 0;
}


