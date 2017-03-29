// This initializes our simulation space
#include <iostream>
#include <vector>
#include "init.hpp"
using namespace std;

void tester()
{
  cout << "foo?" << endl;

}

// Construct and populate the simulation space
// Currently, this creates a symetrical array of atoms
// Obviously, this will need to be improved later
int Init::initialize(int *N, float *sl, float *T, vector<float> *pos, vector<float> *vel)
{
  cout << "Initializer called..." << endl;

  //use this to reference: (*pos).stuff

  tester();
 
  return 0;
}


