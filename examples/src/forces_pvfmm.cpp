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
int Forces_pvfmm::elc_pvfmm(int *N, double *sl, vector<double> *q, vector<double> *pos, vector<double> *force, MPI_Comm* comm, int *start, int *end)
{
  // Variables
  // N		: total number of particles in this task (int)
  // sl		: length of cube (nm, double)
  // q		: 1D vector of charges (int)
  // pos	: 1D vector of positions (nm, double)
  // force	: 1D vector of forces (kg nm/ns^2, double)
  // ke		: electric constant in (kg nm^3/ns^2 , double)
  // comm	: MPI communicator
  // start	: starting index for task (int)
  // end	: ending index for task (int)

  int mult_order=8; //should be taken as an option later

  int i,j,rank,maxrank;
  double ke = 8.89755e18 / pow(*sl,3); 	//in kg nm^3 / ns^2 C^2 sl^3	

  MPI_Comm_rank(*comm,&rank);
  MPI_Allreduce(&rank,&maxrank,1, MPI_INT, MPI_MAX, *comm);  //get maxrank for printing

  //########## setup pvfmm
  const pvfmm::Kernel<double>& kernal_fn=pvfmm::LaplaceKernel<double>::gradient();

  vector<double> dummy(0); //empty vector for the surface

  //create the memory manager (optional?)
  pvfmm::mem::MemoryManager mem_mgr(10000000);

  size_t max_pts=1;
  pvfmm::PtFMM_Tree* tree=PtFMM_CreateTree(*pos, *q, dummy, dummy, *pos, *comm, max_pts, pvfmm::FreeSpace );

  //initialize the matricies, only needs to be done once
  pvfmm::PtFMM matrices(&mem_mgr);
  matrices.Initialize(mult_order, *comm, &kernal_fn);

  //actually build tree
  tree->SetupFMM(&matrices);

  //evaluate tree
  PtFMM_Evaluate(tree, *force, *N); 

  //delete tree once finished
  delete tree;
  //##########

  MPI_Barrier(*comm);
  //Pay no attention to the man behind the curtain!
  for (i=0;i<=maxrank;i++) {
    if (rank == i) {
      ofstream ffile;                 
      if (!rank){
        ffile.open("forces.txt", ofstream::out);
      } else {
        ffile.open("forces.txt", ofstream::out|ofstream::app);
      }
      if (!rank) {
        ffile << "x force, y force, z force (kg nm/ns^2)\n";
      }
      ffile << "output from task : " << rank << "\n";
    for (j=0; j < (*end-*start) ; j++){
      ffile << (*force)[3*j]* *sl << "  " << (*force)[3*j+1]* *sl << "  " << (*force)[3*j+2]* *sl << "\n";
    }
    ffile.close();
  }
  MPI_Barrier(*comm);
}


  return 0;
};

