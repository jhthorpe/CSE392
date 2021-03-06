//header file for forces class
#ifndef _FORCES_PVFMM_HPP_
#define _FORCES_PVFMM_HPP_
#include <vector>
#include <mpi.h>
using namespace std;

class Forces_pvfmm
{
  public:
    //using pvfmm to impliment electrostatics
    //number of mol, box length, charges, position vector, forces vector, mpi communicator, starting index, ending index 
    int elc_pvfmm(int *,double *, vector<double> *, vector<double> *, vector<double> *, MPI_Comm *, int *, int *);
};


#endif
