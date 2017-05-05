// Force evaluations
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <math.h>

//Specific to pvfmm
#include <mpi.h>
#include <pvfmm.hpp>

#include "forces_pvfmm.hpp"
using namespace std;

//~~~     Electrostatic forces. PvFMM implimentation.     ~~~//
int Forces_pvfmm::elc_pvfmm(int *N, double *sl, vector<double> *q, vector<double> *pos, vector<double> *force)
{
  // Variables
  // N		: number of particles (int)
  // sl		: length of cube (nm, int)
  // q		: 1D vector of charges (int)
  // pos	: 1D vector of positions (nm, double)
  // force	: 1D vector of forces (kg nm/ns^2, double)
  // ke		: electric constant in (kg nm^3/ns^2 , double)
  //

  //########## setup pvfmm

  //##########

  
  int i,j;
  double ke = 8.89755e18; 	//in kg nm^3 / ns^2 C^2	

  double fx,fy,fz,f;
  double rx,ry,rz,r;

  for (i=0; i < *N; i++)
  {
    for (j=i+1; j < *N; j++)
    {
      rx=((*pos)[3*j]-(*pos)[3*i]);
      ry=((*pos)[3*j+1]-(*pos)[3*i+1]);
      rz=((*pos)[3*j+2]-(*pos)[3*i+2]);
      r=pow(pow(rx,2)+pow(ry,2)+pow(rz,2),0.5);

      f = ke * (*q)[i]*(*q)[j]/pow(r,2);

      (*force)[3*i] -= f * rx/r;
      (*force)[3*j] += f * rx/r;
      (*force)[3*i+1] -= f * ry/r;
      (*force)[3*j+1] += f * ry/r;
      (*force)[3*i+2] -= f * rz/r;
      (*force)[3*j+2] += f * rz/r;
    };
  };
  //Comment all this out if you don't want to see this 
  ofstream ffile;
  ffile.open("forces.txt");
  ffile << "fx, fy, fz (kg nm / ns^2)\n";
  for (i=0; i < *N ; i++)
  {
    ffile << (*force)[3*i] << "  " << (*force)[3*i+1] << "  " << (*force)[3*i+2] << "\n"; 
  };
  return 0;
};

//~~~     Electrostatic forces. Complete summation. Slow.     ~~~//
