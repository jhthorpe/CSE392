//Header file for energy class
#ifndef _ENERGY_HPP_
#define _ENERGY_HPP_
#include <vector>
using namespace std;

class energy 
{
public:
//Calculate total kinetic energy of system sequenctially: Number of atoms, mass vector, velocity vector, kinetic energy, temperature
  int kinetic_seq(int *, vector<double> *, vector<double> *, double *, double *);

//Calculate Lennard Jones potential energy inside the system sequentially: Number of atoms, side length, sigma, epsilon, position vector, LJ potential energy
  int LJpot_seq(int *, double *, double *, double *, vector<double> *, double *);

//Calculate electrostatic potential energy inside the system sequentially: Number of atoms, charge vector, position vector, electrostatic potential energy
  int elcpot_seq(int *, vector<double> *, vector<double> *, double *);

};

#endif
