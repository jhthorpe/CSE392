// header file for class verlet                                                                                                                                      
#ifndef _VERLET_PAR_HPP_
#define _VERLET_PAR_HPP_
#include <vector>
#include <mpi.h>
using namespace std;

class verlet_par
{
public:
  // number of mol, number of time steps, time step, side length, sigma, epsilon, charge vector, mass vector, position vector, velocity vector, MPI communicator                                                                                                      
  void Integration_mpi(int *, int *, int *, double *, double *, double *, double *, vector<double> *, vector<double> *, vector<double> *, vector<double> *, MPI_Comm *);
};

#endif
