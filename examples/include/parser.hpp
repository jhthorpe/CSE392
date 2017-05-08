// header file for our Parser class
#ifndef _PARSER_HPP_
#define _PARSER_HPP_
#include <mpi.h>

class Parser
{
  public:
    //num part, box length, temp, mass, timestep, num steps, LJ sigma, LJ epsilon, charge, MPI compiler
  int getInput(int *, double * ,double *, double *, double *, int*, int *, double *, double *, double *, MPI_Comm * );	//actual parser
};
#endif
