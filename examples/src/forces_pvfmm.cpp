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
  // sl		: length of cube (nm, int)
  // q		: 1D vector of charges (int)
  // pos	: 1D vector of positions (nm, double)
  // force	: 1D vector of forces (kg nm/ns^2, double)
  // ke		: electric constant in (kg nm^3/ns^2 , double)
  //

  int mult_order=2; //should be taken as an option later

  int i,j;
  double ke = 8.89755e18; 	//in kg nm^3 / ns^2 C^2	

  double fx,fy,fz,f;
  double rx,ry,rz,r;

  cout << "here1" << endl;

  //########## setup pvfmm
  const pvfmm::Kernel<double>& kernal_fn=pvfmm::LaplaceKernel<double>::gradient();

  vector<double> dummy(0); //empty vector for the surface

  //change distances (here for now)
  for (i=0;i < *N;i++){
    (*pos)[3*i] = (*pos)[3*i] / *sl;
    (*pos)[3*i+1] = (*pos)[3*i+1] / *sl;
    (*pos)[3*i+2] = (*pos)[3*i+2] / *sl;
  } 

  pvfmm::mem::MemoryManager mem_mgr(10000000);

  cout << "before tree" << endl;

  std::cout<<pos->size()<<' '<<q->size()<<'\n';
  
  size_t max_pts=1;
  pvfmm::PtFMM_Tree* tree=PtFMM_CreateTree(*pos, *q, dummy, dummy, *pos, *comm, max_pts, pvfmm::FreeSpace );

  cout << "before eval" << endl;

  pvfmm::PtFMM matrices(&mem_mgr);
  matrices.Initialize(mult_order, *comm, &kernal_fn);

  cout << "before eval" << endl;

  tree->SetupFMM(&matrices);

  cout << "before eval" << endl;

  PtFMM_Evaluate(tree, *force, *N); 

  cout << "after eval" << endl;
  delete tree;
  //##########

  

  //Comment all this out if you don't want to see this 
  ofstream ffile;
  ffile.open("forces.txt");
  ffile << "fx, fy, fz (kg nm / ns^2)\n";
  for (i=0; i < *N ; i++)
  {
    ffile << (*force)[3*i] << "  " << (*force)[3*i+1] << "  " << (*force)[3*i+2] << "\n"; 
  };

  return 0;
};

//~~~     Electrostatic forces. Complete summation. Slow.     ~~~//
