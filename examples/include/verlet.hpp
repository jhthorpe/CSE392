// header file for class verlet
                                                                            
#ifndef _VERLET_HPP_
#define _VERLET_HPP_
#include <mpi.h>
#include <vector>
using namespace std;

class verlet
{
public:
  // number of mol, number of time steps, time step, side length, sigma, epsilon, charge vector, mass vector, position vector, velocity vector, MPI Communicator, start index, end index
                                                                                     
  void Integration(int *, int *, int *, double *, double *, double *, double *, vector<double> *, vector<double> *, vector<double> *, vector<double> *, MPI_Comm *, int *, int *);
};

#endif

