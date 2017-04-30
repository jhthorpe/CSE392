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
  float kB = 1.38064852e-23; 		// in J/K

  cout << "Initializer called..." << endl;

  int pl = ceil(pow(*N,1.0/3.0));
  float bl = float(*sl /float(pl+1)); 
  float stdv = pow(kB * *T / *m,0.5);

  random_device generator;
  normal_distribution<float> vdist(0,stdv);
 
  // Create our cube of molecules 
  for (i=0; i < *N; i++) 
  {
    // I could, of course, have used 3 loops, but this is better flops/mops ;)

    (*pos)[i] = bl * float(i % pl);						// x position
    (*pos)[i+1] = bl * float((i / pl) % pl);					// y position
    (*pos)[i+2] = bl * float((i / int(pow(pl,2.0))) % int(pow(pl,2.0)));	// z position 

    // st. dev = sqrt(kB * T / m)
    // call std::normal_distribution(0,st.dev), as average velocity is 0
    

    //create random number generator and distribution object

    (*vel)[i]=vdist(generator);
    (*vel)[i+1]=vdist(generator);				// y velocity
    (*vel)[i+2]=vdist(generator);				// z velocity

    cout << i << endl;
    cout << (*pos)[i] << ", " << (*pos)[i+1] << ", " << (*pos)[i+2] << endl;
    cout << (*vel)[i] << ", " << (*vel)[i+1] << ", " << (*vel)[i+2] << endl;
  }
 
  return 0;
}


