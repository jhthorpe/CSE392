// header file for class verlet                                                                                                                                      
#ifndef _VERLET_HPP_
#define _VERLET_HPP_
#include <vector>
using namespace std;

class verlet
{
public:
  // number of mol, number of time steps, time step, side length, sigma, epsilon, charge vector, mass vector, position vector, velocity vector                                                                                                       
  void Integration(int *, int *, int *, double *, double *, double *, double *, vector<double> *, vector<double> *, vector<double> *, vector<double> *);
};

#endif