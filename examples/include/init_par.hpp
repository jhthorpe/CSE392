// header file for our Init (initializer) class
#ifndef _INIT_PAR_HPP_
#define _INIT_PAR_HPP_
#include <vector>
#include <mpi.h>
using namespace std;

class Init_par
{
  public:
    // number of mol, box length, temp, mass, position vector, velocity vector, starting index, ending index
    int initialize_mpi(int *, double *,double *, vector<double> *, vector<double> *, vector<double> * , int *, int *, MPI_Comm *);
};

#endif

