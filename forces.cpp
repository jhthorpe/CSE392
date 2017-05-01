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

//sequential LJ force evaluation, no boundary conditions. Slow.
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
  int i,j;
  float f, r,rx,ry,rz;
  cout << "Leonard Jones called" << endl;
 // r=1.0;
 // f = (24e0 * *eps / *sig) * (2e0 * pow(*sig/r,13) - pow(*sig/r,7) );
 // cout << "single force over r nm : " << f << endl;

  //check that multiplying by unit vector is alright 
  for (i=0; i < *N; i++)
  {
    for (j=i+1; j < *N; j++)
    { 
      rx=((*pos)[3*j]-(*pos)[3*i]);
      ry=((*pos)[3*j+1]-(*pos)[3*i+1]);
      rz=((*pos)[3*j+2]-(*pos)[3*i+2]);
      r=pow(pow(rx,2)+pow(ry,2)+pow(rz,2),0.5);

      f = (24e0 * *eps / *sig) * (2e0 * pow(*sig/r,13) - pow(*sig/r,7) );
//      cout << "force i= " << i << ", j= "<< j << " : " << f << endl;

     (*force)[3*i] += f * (rx/r);
     (*force)[3*j] -= f * (rx/r);
     (*force)[3*i+1] += f * (ry/r);
     (*force)[3*j+1] -= f * (ry/r);
     (*force)[3*i+2] += f * (rz/r);
     (*force)[3*j+2] -= f * (rz/r);
 
//      (*force)[3*i] += (24e0 * *eps / *sig) * (2e0 * pow(*sig/rx,13) - pow(*sig/rx,7) )*(rx/r);	// x force in kg nm/ ns^2	
//      (*force)[3*j] -= (24e0 * *eps / *sig) * (2e0 * pow(*sig/rx,13) - pow(*sig/rx,7) )*(rx/r);	// x force in kg nm/ ns^2	
//      (*force)[3*i+1] += (24e0 * *eps / *sig) * (2e0 * pow(*sig/ry,13) - pow(*sig/ry,7) )*(ry/r);	// y force in kg nm/ ns^2	
//      (*force)[3*j+1] -= (24e0 * *eps / *sig) * (2e0 * pow(*sig/ry,13) - pow(*sig/ry,7) )*(ry/r);	// y force in kg nm/ ns^2	
//      (*force)[3*i+2] += (24e0 * *eps / *sig) * (2e0 * pow(*sig/rz,13) - pow(*sig/rz,7) )*(rz/r);	// z force in kg nm/ ns^2	
//      (*force)[3*j+2] -= (24e0 * *eps / *sig) * (2e0 * pow(*sig/rz,13) - pow(*sig/rz,7) )*(rz/r);	// z force in kg nm/ ns^2	
    };
  };
 
  cout << "forces on particle 0:" << (*force)[0] <<","<<(*force)[1]<<","<<(*force)[2]<<","<< endl; 
  cout << "forces on particle 1:" << (*force)[3] <<","<<(*force)[4]<<","<<(*force)[5]<<","<< endl; 
  cout << "forces on particle 2:" << (*force)[6] <<","<<(*force)[7]<<","<<(*force)[8]<<","<< endl;


  return 0;  
};
