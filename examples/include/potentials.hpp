//header file for forces class
#ifndef _POTENTIALS_HPP_
#define _POTENTIALS_HPP_
#include <vector>
#include <mpi.h>
using namespace std;

class potentials
{
public:
  //using pvfmm to evaluate electrostatics and lennard jones potential and forces

  //number of mol, box length, charges, position vector, forces vector, electric potential energy, mpi communicator, task start index
  int elc_pvfmm(int *,double *, vector<double> *, vector<double> *, vector<double> *, double *, MPI_Comm *, int *);

  //number of mol, box length, charges, position vector, forces vector, LJ potential energy
  int LJ_pvfmm(int *, double *, double *, double *, vector<double> *, vector<double> *, double *, MPI_Comm *);
};


#endif

