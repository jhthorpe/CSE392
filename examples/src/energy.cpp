// Energy evaluations                                                                                                                                                           
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <math.h>
#include "energy.hpp"
using namespace std;


//Evaluate kinetic energy of all atoms

int energy::kinetic_seq(int *N, vector<double> *mass, vector<double> *vel, double *kin, double *temp)
{
  //Variables
  // N          : number of input molecules
  // sl         : side length of box (nm)
  // mass       : 1D vector of masses of particles (size N)                                                                                                                   
  // vel        : 1D vector of particle velocities (size 3N)
  // kin        : Total kinetic energy (kg nm^2/ns^2)
  // temp       : Temperature in (degrees Kelvin)

  double kB = 1.38064852e-23;     //Boltzmann constant (kg-nm^2/ns^2-K)

  //internal variables                                                                                                                                                  
  int i;

  for (i=0; i < *N; i++)
    {
      *kin += 0.5 * ((*mass)[i]) * ( pow((*vel)[3*i],2) + pow((*vel)[3*i+1],2) + pow((*vel)[3*i+2],2) ) ;
    }

  *temp = (2.0/3.0) * (*kin)/(*N*kB);

  return 0;
}

int energy::LJpot_seq(int *N, double *sl, double *sig, double *eps, vector<double> *pos, double *LJpot)
{
  //Variables
  // N          : number of input molecules
  // sl         : side length of box (nm)
  // sig        : sigma in LJ potential
  // eps        : epsilon in LJ potential
  // pos        : 1D vector of particle positions (size 3N)
  // LJpot      : Total LJ potential energy

  //internal variables
  int i,j;
  double r,rx,ry,rz;

  for (i=0; i < *N; i++)
    {
      for (j=i+1; j < *N; j++)
	{
	  rx=((*pos)[3*j]-(*pos)[3*i]);
	  ry=((*pos)[3*j+1]-(*pos)[3*i+1]);
	  rz=((*pos)[3*j+2]-(*pos)[3*i+2]);
	  r=pow(pow(rx,2)+pow(ry,2)+pow(rz,2),0.5);

	  *LJpot += (4e0 * *eps / *sig) * (pow(*sig/r,12) - pow(*sig/r,6) );
	};
    };
  return 0;
}

int energy::elcpot_seq(int *N, vector<double> *q, vector<double> *pos, double *elcpot) 
{
  //Variables
  // N          : number of atoms
  // q          : 1D vector of charges (size N)
  // pos        : 1D vector of particle positions (size 3N)
  // elcpot     : Total electrostatic  potential energy

  //Permittivity
  double ke = 8.89755e18;       //in kg nm^3 / ns^2 C^2 

  //internal variables
  int i,j;
  double r,rx,ry,rz;

  for (i=0; i < *N; i++)
    {
      for (j=i+1; j < *N; j++)
        {
          rx=((*pos)[3*j]-(*pos)[3*i]);
          ry=((*pos)[3*j+1]-(*pos)[3*i+1]);
          rz=((*pos)[3*j+2]-(*pos)[3*i+2]);
          r=pow(pow(rx,2)+pow(ry,2)+pow(rz,2),0.5);

	  *elcpot += ke * (*q)[i]*(*q)[j]/r;
        };
    };
  return 0;
}