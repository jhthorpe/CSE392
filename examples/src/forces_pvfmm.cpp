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
int Forces_pvfmm::elc_pvfmm(int *N, double *sl, vector<double> *q, vector<double> *pos, vector<double> *force, MPI_Comm* comm)
{
  // Variables
  // N		: number of particles (int)
  // sl		: length of cube (nm, double)
  // q		: 1D vector of charges (int)
  // pos	: 1D vector of positions (nm, double)
  // force	: 1D vector of forces (kg nm/ns^2, double)
  // ke		: electric constant in (kg nm^3/ns^2 , double)
  //

  int mult_order=2; //should be taken as an option later

  int i;
  double ke = 8.89755e18 / pow(*sl,3); 	//in kg nm^3 / ns^2 C^2 sl^3	

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

  int rank;
  MPI_Comm_rank(*comm,&rank);
  if (rank == 0 || rank == 1) {
    ofstream ffile;
    if (!rank) {
      //ofstream ffile0;
      ffile.open("task0");
    } else {
      //ofstream ffile1;
      ffile.open("task1");
    }
    ffile << "forces from task " << rank << endl;
    for (i=0; i < *N ; i++) ffile << (*force)[3*i]* *sl << "  " << (*force)[3*i+1]* *sl << "  " << (*force)[3*i+2]* *sl << "\n"; 
    ffile.close();
  }

  //Comment all this out if you don't want to see this 
  /*
  ofstream ffile;
  ffile.open("forces.txt");
  ffile << "fx, fy, fz (kg nm / ns^2)\n";
  for (i=0; i < *N ; i++)
  {
    ffile << (*force)[3*i]* *sl << "  " << (*force)[3*i+1]* *sl << "  " << (*force)[3*i+2]* *sl << "\n"; 
  };
  ffile.close();
  */

  return 0;
};

