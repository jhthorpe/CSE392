// Force and Energy evaluations
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <math.h>
#include <mpi.h>

#include <pvfmm.hpp>
#include "potentials.hpp"
using namespace std;

//Global variables, to be assigned values upon function calls
double eps_global, sig_global, ke_global;

//Laplace potential kernel
template <class Real>
void laplace_poten(Real* r_src, int src_cnt, Real* v_src, int dof, Real* r_trg, int trg_cnt, Real* v_trg, pvfmm::mem::MemoryManager* mem_mgr) {
  const Real OOFP = ke_global;
  for (int t = 0; t < trg_cnt; t++) {
    Real p = 0;
    for (int s = 0; s < src_cnt; s++) {
      Real dR[3] = {r_trg[3 * t    ] - r_src[3 * s    ],
                    r_trg[3 * t + 1] - r_src[3 * s + 1],
                    r_trg[3 * t + 2] - r_src[3 * s + 2]};
      Real R2 = (dR[0] * dR[0] + dR[1] * dR[1] + dR[2] * dR[2]);
      if (R2 != 0) p += v_src[s] / sqrt(R2);
    }
    v_trg[t] += p * OOFP;
  }
}

//Laplace potential and gradient function for electrostatic interaction. This will be used to build kernel for evaluations using PVFMM
template <class Real>
void laplace_poten_and_grad(Real* r_src, int src_cnt, Real* v_src, int dof, Real* r_trg, int trg_cnt, Real* v_trg, pvfmm::mem::MemoryManager* mem_mgr) {
  const Real OOFP = ke_global;
  for (int t = 0; t < trg_cnt; t++) {
    Real p[4] = {0.0,0.0,0.0,0.0};
    for (int s = 0; s < src_cnt; s++) {
      Real dR[3] = {r_trg[3 * t    ] - r_src[3 * s    ],
                    r_trg[3 * t + 1] - r_src[3 * s + 1],
                    r_trg[3 * t + 2] - r_src[3 * s + 2]};
      Real R2 = (dR[0] * dR[0] + dR[1] * dR[1] + dR[2] * dR[2]);   //Distance squared between r_trg and r_src
      if (R2 != 0){
         Real invR = 1.0/sqrt(R2);                                  //Inverse distance
         Real invR3=invR*invR*invR;                                 //Inverse Distance cubed
         p[0] += v_src[s] * invR3 * dR[0];                          //Force in x-direction
         p[1] += v_src[s] * invR3 * dR[1];                          //Force in y-direction
         p[2] += v_src[s] * invR3 * dR[2];                          //Force in z-direction
         p[3] += v_src[s] * invR;                                   //Potential energy
       }
     }
     for (int i=0;i<4;i++) {
       v_trg[t*4+i] += p[i] * OOFP;                                 //Multiply by electric constant to get correct units
       //       cout << "elc v_trg[t*4+" << i << "] =" << v_trg[t*4+i] << endl;
     }
   }
 }
 
 //Lennard Jones potential and gradient function. This will be used to build kernel for evaluations using PVFMM
 template <class Real>
 void lj_poten_and_grad(Real* r_src, int src_cnt, Real* v_src, int dof, Real* r_trg, int trg_cnt, Real* v_trg, pvfmm::mem::MemoryManager* mem_mgr) {
   const Real LJFP = 4.0 * eps_global;
   for (int t = 0; t < trg_cnt; t++) {
     Real p[4] = {0.0,0.0,0.0,0.0}, q[4] = {0.0,0.0,0.0,0.0};
     for (int s = 0; s < src_cnt; s++) {
       Real dR[3] = {r_trg[3 * t    ] - r_src[3 * s    ],
                     r_trg[3 * t + 1] - r_src[3 * s + 1],
                     r_trg[3 * t + 2] - r_src[3 * s + 2]};
       Real R2 = (dR[0] * dR[0] + dR[1] * dR[1] + dR[2] * dR[2]);   //Distance squared between r_trg and r_src
       if (R2 != 0){
         Real invR2 = 1.0/R2;                                       //Inverse distance squared (scaled by sigma)
         Real invR12= pow(invR2,6), sig12= pow(sig_global,12); 
	 Real invR6 = pow(invR2,3), sig6 = pow(sig_global,6 );  
         p[0] += 12 * v_src[s] * sig12 * invR12 * invR2 * dR[0];
         q[0] += 6  * v_src[s] * sig6  *  invR6 * invR2 * dR[0];            //Forces in x-direction
         p[1] += 12 * v_src[s] * sig12 * invR12 * invR2 * dR[1];
         q[1] += 6  * v_src[s] * sig6  *  invR6 * invR2 * dR[1];            //Forces in y-direction
         p[2] += 12 * v_src[s] * sig12 * invR12 * invR2 * dR[2];
         q[2] += 6  * v_src[s] * sig6  *  invR6 * invR2 * dR[2];            //Forces in z-direction
         p[3] +=      v_src[s] * sig12 * invR12;
         q[3] +=      v_src[s] * sig6  *  invR6;                            //Potential energies
       }
     }
     for (int i=0;i<4;i++) {
       v_trg[t*4+i] += (p[i] - q[i]) * LJFP;                        //Final values of LJ forces and energies
       //       cout << "LJ v_trg[t*4+" << i << "] =" << v_trg[t*4+i] << endl;
     }
   }
}

//~~~     Electrostatic forces and potential energy. PvFMM implementation.     ~~~//
int potentials::elc_pvfmm(int *N, double *sl, vector<double> *q, vector<double> *pos, vector<double> *force, double *elcpoten_energy, MPI_Comm* comm, int *)
{
  cout << "~~~~~~~~~~~~~~elc_pvfmm called~~~~~~~~~~~~~" << endl;
  // Variables                                                                                                                               
  // N               : number of particles (int)
  // sl              : size length of box (nm, double)
  // q               : 1D vector of charges (C, double)
  // pos             : 1D vector of positions (nm / sl, double)
  // force           : 1D vector of forces (kg nm/ns^2 / sl, double)
  // elcpoten_energy : electrostatic potential energy of system (kg nm^2/ns^2 / sl^2)
  // comm            : pvfmm related MPI communicator       
  // poten_and_grad  : 1D vector of forces and potential (double)
  // ke              : electric constant in (ideally in kg nm^3/ns^2 , double)
  // comm	     : MPI communicator
  // *	     : *ing index of this task

  int mult_order=8;                                         //8 or 10 for electrostatic interactions
  int i;
  double ke = 8.89755e18 * (*q)[0] * (*q)[0] / pow(*sl,3);  //in kg nm^3 / ns^2 sl^3; scaled with respect to charges for kernel defintion. THIS REQUIRES ALL CHARGES TO BE SAME 
  ke_global = ke;
  vector<double> src_value(*N, 1.0);                        //default value for src_value 
  vector<double> poten_and_grad(4 * *N, 0.0);               //vector to store output from laplace_poten_and_grad()
  vector<double> dummy(0);                                  //empty vector for the surface

  //Build kernels                                                                                                       
  const pvfmm::Kernel<double> potn_ker=pvfmm::BuildKernel<double, laplace_poten<double> >("laplace", 3, std::pair<int,int>(1,1), NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL);
  const pvfmm::Kernel<double> elc_ker=pvfmm::BuildKernel<double, laplace_poten_and_grad<double> >("laplace-poten-and-grad", 3, std::pair<int,int>(1,4), &potn_ker ,&potn_ker,NULL,&potn_ker ,&potn_ker ,NULL ,&potn_ker ,NULL);

  //create the memory manager
  pvfmm::mem::MemoryManager mem_mgr(10000000);

  size_t max_pts=1;
  pvfmm::PtFMM_Tree* tree=PtFMM_CreateTree(*pos, src_value, dummy, dummy, *pos, *comm, max_pts, pvfmm::Periodic );

  //initialize the matrices, only needs to be done once
  pvfmm::PtFMM matrices(&mem_mgr);
  matrices.Initialize(mult_order, *comm, &elc_ker);

  //actually build tree
  tree->SetupFMM(&matrices);

  //evaluate tree
  pvfmm::PtFMM_Evaluate(tree, poten_and_grad, *N);

  cout << "Values of first 4 numbers in poten_and_grad are: " << poten_and_grad[0] << ", " << poten_and_grad[1] << ", " << poten_and_grad[2] << " & " << poten_and_grad[3] << endl;

  //delete tree once finished
  delete tree;

  //Print output
  int rank;
  MPI_Comm_rank(*comm, &rank);
  int maxrank,np;
  MPI_Comm_size(*comm, &np);
  MPI_Allreduce(&rank,&maxrank,1, MPI_INT, MPI_MAX, *comm);


  for (int j=0; j < *N ; j++){
    (*force)[3*(j) ]  += poten_and_grad[4*(j)  ];
    (*force)[3*(j)+1]  += poten_and_grad[4*(j)+1];
    (*force)[3*(j)+2]  += poten_and_grad[4*(j)+2];
    *elcpoten_energy += poten_and_grad[4*(j)+3];
  }

  MPI_Barrier(*comm);

  for (i=0;i<=maxrank;i++){
    if (rank==i){
      ofstream ffile;
      if (!rank) {
        ffile.open("forces.txt", ofstream::out);
        ffile << "fx, fy, fz (kg nm / ns^2)\n";
      } else {
        ffile.open("forces.txt", ofstream::out|ofstream::app);
      }
      ffile << "hi from task " << rank << "\n";
      for (int j=0; j < *N ; j++){
        ffile << (*force)[3*(j)]* *sl << "  " << (*force)[3*(j)+1]* *sl << "  " << (*force)[3*(j)+2]* *sl << "\n"; //All forces are re-scaled back to Newtons by muliplying by side length
      }
      ffile.close();
    }
    MPI_Barrier(*comm);
  }
  return 0;
};

//~~~     Lennard-Jones forces and potential energy. PvFMM implementation.     ~~~//
int potentials::LJ_pvfmm(int *N, double *sl, double *sig, double *eps, vector<double> *pos, vector<double> *force, double *LJpoten_energy, MPI_Comm* comm)
{
  cout << "~~~~~~~~~~~~~~LJ_pvfmm called~~~~~~~~~~~~~" << endl;
  // Variables
  // N             : number of particles (int)
  // sl            : size length of box (nm)
  // sig           : sigma in LJ potential (nm / sl, double)
  // eps           : epsilon in LJ potential (kg-nm^2/ns^2 / sl^2, double)
  // pos           : 1D vector of positions (nm / sl, double)
  // force         : 1D vector of forces (kg nm/ns^2 / sl, double)
  // LJpoten_energy: Lennard-Jones potential energy of system (kg nm^2/ns^2 / sl^2)
  // poten_and_grad: 1D vector of forces and potential (double)
  // comm          : pvfmm related MPI communicator
  // src_value     : pvfmm related source value for each pairwise interaction (all equal to 1.0)

  int mult_order=0;                             //0 for lennard jones
  int i;
  vector<double> src_value(*N, 1.0);            //default value for src_value
  vector<double> poten_and_grad(4 * *N, 0.0);   //Output vector for lj_poten_and_grad()

  eps_global = *eps; sig_global = *sig;             

  //Build kernels
  const pvfmm::Kernel<double> lj_ker=pvfmm::BuildKernel<double, lj_poten_and_grad<double> >("lj-poten-and-grad", 3, std::pair<int,int>(1,4));

  vector<double> dummy(0); //empty vector for the surface

  //create the memory manager
  pvfmm::mem::MemoryManager mem_mgr(10000000);

  size_t max_pts = 10;
  pvfmm::PtFMM_Tree* tree=PtFMM_CreateTree( *pos, src_value, dummy, dummy, *pos, *comm, max_pts, pvfmm::Periodic );

  //initialize the matricies, only needs to be done once
  pvfmm::PtFMM matrices(&mem_mgr);
  matrices.Initialize(mult_order, *comm, &lj_ker);

  //actually build tree
  tree->SetupFMM(&matrices);

  //evaluate tree
  PtFMM_Evaluate(tree, poten_and_grad, *N);

  cout << "Values of first 4 numbers in poten_and_grad are: " << poten_and_grad[0] << ", " << poten_and_grad[1] << ", " << poten_and_grad[2] << " & " << poten_and_grad[3] << endl;

  //delete tree once finished
  delete tree;

  //Comment all this out if you don't want to see this
  
  ofstream ljfile;
  ljfile.open("LJforces.txt");
  ljfile << "fx, fy, fz (kg nm / ns^2)\n";
  for (i=0; i < *N ; i++)
    {
      (*force)[3*i  ] += poten_and_grad[4*i  ];
      (*force)[3*i+1] += poten_and_grad[4*i+1];
      (*force)[3*i+2] += poten_and_grad[4*i+2];
      *LJpoten_energy += poten_and_grad[4*i+3];                                         //Final value of potential energy has a correction factor due to kernel definition
      ljfile << (*force)[3*i]* *sl << "  " << (*force)[3*i+1]* *sl << "  " << (*force)[3*i+2]* *sl << "\n"; //All forces re-scaled back to Newton by multiplying by side length 
    };
  ljfile.close();

  return 0;
};
