// Force evaluations
#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include "forces.hpp"
using namespace std;

//tester function
void Forces::Test()
{ 
  cout << "forces called" << endl;
}

//sequential LJ force evaluation
int Forces::LJ_seq(int *N, float *sl, float *sig, float *eps, vector<float> *pos, vector<float> *force)
{
  //Varaibles
  // N		: number of input molecules
  // sl		: side length of box (nm)
  // sig	: LJ sigma in (nm)
  // eps	: LJ epsilon in (kg nm^2/ns^2)
  // pos	: 1D vector of particle positions (in nm)
  // force	: 1D vector of particle forces (in kg nm / ns^2)

  //internal variables
  float f, r;
  r = 0.05; 
  cout << "Leonard Jones called" << endl;

  f = (24e0 * *eps / *sig) * (2e0 * pow(*sig/r,13) - pow(*sig/r,7) );
  cout << "single force over r nm : " << f << endl;
  

  return 0;  
};
