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
double eps_global, ke_global;

//Laplace potential and gradient function for electrostatic interaction. This will be used to build kernel for evaluations using PVFMM
template <class Real>
void laplace_poten_and_grad(Real* r_src, int src_cnt, Real* v_src, int dof, Real* r_trg, int trg_cnt, Real* v_trg, pvfmm::mem::MemoryManager* mem_mgr) {
  Real OOFP = ke_global;
  //  cout << "v_src[0]: " << v_src[0] << endl; 
  for (int t = 0; t < trg_cnt; t++) {
    Real p[4] = {0.0,0.0,0.0,0.0};
    for (int s = 0; s < src_cnt; s++) {
      Real dR[3] = {r_trg[3 * t    ] - r_src[3 * s    ],
                    r_trg[3 * t + 1] - r_src[3 * s + 1],
                    r_trg[3 * t + 2] - r_src[3 * s + 2]};
      Real R2 = (dR[0] * dR[0] + dR[1] * dR[1] + dR[2] * dR[2]);   //Distance squared between r_trg and r_src
      //      cout << "~~~~~~~~~~~~~~laplace_poten_and_grad() called~~~~~~~~~~~~~" << endl;
      //      cout << "Value of R2 = " << R2 << endl; 
      if (R2 != 0){
	//cout << "Value of r_src[s]: " << r_src[s] << endl;
	//cout << "Value of r_trg[s]: " << r_trg[s] << endl;
        Real invR = 1.0/sqrt(R2);                                  //Inverse distance
        Real invR3=invR*invR*invR;                                 //Inverse Distance cubed
        p[0] += v_src[s] * invR3 * dR[0];               //Force in x-direction
        p[1] += v_src[s] * invR3 * dR[1];               //Force in y-direction
        p[2] += v_src[s] * invR3 * dR[2];               //Force in z-direction
        p[3] += v_src[s] * invR;                        //Potential energy
      }
    }
    for (int i=0;i<4;i++) {
      v_trg[t*4+i] += p[i] * OOFP;             //Multiply by electric constant to get correct units
      //      cout << "v_trg[t*4+" << i << "] =" << v_trg[t*4+i] << endl;
    }
  }
}

//Lennard Jones potential and gradient function. This will be used to build kernel for evaluations using PVFMM
template <class Real>
void LJ_poten_and_grad(Real* r_src, int src_cnt, Real* v_src, int dof, Real* r_trg, int trg_cnt, Real* v_trg, pvfmm::mem::MemoryManager* mem_mgr) {
  Real OOFP = 4.0 * eps_global;
  for (int t = 0; t < trg_cnt; t++) {
    Real p[4] = {0.0,0.0,0.0,0.0}, q[4] = {0.0,0.0,0.0,0.0};
    for (int s = 0; s < src_cnt; s++) {
      Real dR[3] = {r_trg[3 * t    ] - r_src[3 * s    ],
                    r_trg[3 * t + 1] - r_src[3 * s + 1],
                    r_trg[3 * t + 2] - r_src[3 * s + 2]};
      Real R2 = (dR[0] * dR[0] + dR[1] * dR[1] + dR[2] * dR[2]);       //Distance squared between r_trg and r_src
      //      cout << "~~~~~~~~~~~~~~LJ_poten_and_grad() called~~~~~~~~~~~~~" << endl;
      //      cout << "Value of R2 = " << R2 << endl; 
      if (R2 != 0){
        Real invR2 = 1.0/R2;                              //Inverse distance squared (scaled by sigma)
        Real invR12= pow(invR2,6), invR6= pow(invR2,3);   //Inverse scaled distance to powers 12 and 6
        p[0] += 12 * v_src[s] * invR12 * invR2 * dR[0]; 
	q[0] += 6  * v_src[s] *  invR6 * invR2 * dR[0];   //Forces in x-direction
        p[1] += 12 * v_src[s] * invR12 * invR2 * dR[1]; 
	q[1] += 6  * v_src[s] *  invR6 * invR2 * dR[1];   //Forces in y-direction
        p[2] += 12 * v_src[s] * invR12 * invR2 * dR[2]; 
	q[2] += 6  * v_src[s] *  invR6 * invR2 * dR[2];   //Forces in z-direction
        p[3] +=      v_src[s] * invR12; 
	q[3] +=      v_src[s] *  invR6;                   //Potential energies
      }
    }
    for (int i=0;i<4;i++) {
      v_trg[t*4+i] += (p[i] - q[i]) * OOFP;        //Final values of LJ forces and energies
      //      cout << "v_trg[t*4+" << i << "] =" << v_trg[t*4+i] << endl;
    }
  }
}

//Build kernels for each of these functions
const pvfmm::Kernel<double> elc_ker=pvfmm::BuildKernel<double, laplace_poten_and_grad<double> >("laplace-poten-and-grad", 3, std::pair<int,int>(1,4), &elc_ker ,&elc_ker ,NULL,&elc_ker ,&elc_ker ,NULL ,&elc_ker ,NULL);
const pvfmm::Kernel<double> LJ_ker=pvfmm::BuildKernel<double, LJ_poten_and_grad<double> >("LJ-poten-and-grad", 3, std::pair<int,int>(1,4), &LJ_ker ,&LJ_ker ,NULL,&LJ_ker ,&LJ_ker ,NULL ,&LJ_ker ,NULL);

//~~~     Electrostatic forces and potential energy. PvFMM implementation.     ~~~//
int potentials::elc_pvfmm(int *N, double *sl, vector<double> *q, vector<double> *pos, vector<double> *force, double *elcpoten_energy, MPI_Comm* comm)
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

  int mult_order=10; //8 or 10 for electrostatic interactions
  int i;
  double ke = 8.89755e18 * (*q)[0] * (*q)[0] / pow(*sl,3);  //in kg nm^3 / ns^2 sl^3; scaled with respect to charges for kernel defintion. THIS REQUIRES ALL CHARGES TO BE SAME 
  ke_global = ke;
  vector<double> poten_and_grad(4 * *N, 0.0);

  //########## setup pvfmm                              
  //  const pvfmm::Kernel<double>& kernal_fn=elc_ker<double>::laplace_poten_and_grad();
  //  const pvfmm::Kernel<double>& kernal_fn=elc_ker;
  //  const pvfmm::Kernel<double>& kernal_fn=laplace_poten_and_grad(*pos,N,*q,1,*pos,&poten_and_grad,comm);
  vector<double> dummy(0); //empty vector for the surface

  //create the memory manager (optional?)
  pvfmm::mem::MemoryManager mem_mgr(10000000);

  size_t max_pts=1;
  pvfmm::PtFMM_Tree* tree=PtFMM_CreateTree(*pos, *q, dummy, dummy, *pos, *comm, max_pts, pvfmm::FreeSpace );

  //initialize the matricies, only needs to be done once
  pvfmm::PtFMM matrices(&mem_mgr);
  matrices.Initialize(mult_order, *comm, &elc_ker);

  //actually build tree
  tree->SetupFMM(&matrices);

  //evaluate tree
  PtFMM_Evaluate(tree, poten_and_grad, *N);
  cout << "~~~~~~~~~~~~~~elc_pvfmm tree evaluated~~~~~~~~~~~~~" << endl;
  cout << "Values of first 4 numbers in poten_and_grad are: " << poten_and_grad[0] << ", " << poten_and_grad[1] << ", " << poten_and_grad[2] << " & " << poten_and_grad[3] << endl;
  //delete tree once finished
  delete tree;
  
  //##########

  //Comment all this out if you don't want to see this
  ofstream ffile;
  ffile.open("forces.txt"); 
  ffile << "fx, fy, fz (kg nm / ns^2)\n";
  for (i=0; i < *N ; i++)
    {
      (*force)[3*i  ]  += poten_and_grad[4*i  ];
      (*force)[3*i+1]  += poten_and_grad[4*i+1];
      (*force)[3*i+2]  += poten_and_grad[4*i+2];
      *elcpoten_energy += poten_and_grad[4*i+3];
      ffile << (*force)[3*i]* *sl << "  " << (*force)[3*i+1]* *sl << "  " << (*force)[3*i+2]* *sl << "\n"; //All forces are re-scaled back to units of (N) by multiplying by side length
    };

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

  int mult_order=0; //0 for lennard jones
  int i;
  vector<double> src_value(3 * *N, 1.0); 
  //  cout << "src_values: " << src_value[0] << ", " << src_value[1] << ", " << src_value[2] << ", " << src_value[3] << ", " << src_value[4] << ", " << src_value[5] << endl;
  eps_global = *eps / *sig ;                    //For the sake of kernel function
  cout << "eps_global: " << eps_global << endl;
  vector<double> poten_and_grad(4 * *N, 0.0);   //No unique units
  vector<double> scaled_pos(3 * *N);            //Final units (nm / sl) * (sl / nm) 
  for (i=0; i< (3 * *N); i++) {
    scaled_pos[i] = (*pos)[i] / *sig ;       //Scale position with respect to sigma for LJ kernel calculation
    cout << "pos[" << i << "]: "<< (*pos)[i] << endl;
    cout << "scaled_pos[" << i << "]: " << scaled_pos[i] << endl;
  }
  //  cout << "pos: " << (*pos)[0] << ", " << (*pos)[1] << ", " << (*pos)[2] << ", " << (*pos)[3] << ", " << (*pos)[4] << ", " << (*pos)[5] << endl;
  //  cout << "scaled_pos: " << scaled_pos[0] << ", " << scaled_pos[1] << ", " << scaled_pos[2] << ", " << scaled_pos[3] << ", " << scaled_pos[4] << ", " << scaled_pos[5] << endl;

  //########## setup pvfmm
  //  const pvfmm::Kernel<double>& kernal_fn=LJ_ker<double>::LJ_poten_and_grad();
  //  const pvfmm::Kernel<double>& kernal_fn=LJ_ker;
  vector<double> dummy(0); //empty vector for the surface

  //create the memory manager (optional?)
  pvfmm::mem::MemoryManager mem_mgr(10000000);

  size_t max_pts=1;
  pvfmm::PtFMM_Tree* tree=PtFMM_CreateTree( scaled_pos, src_value, dummy, dummy, scaled_pos, *comm, max_pts, pvfmm::FreeSpace );
  cout << "~~~~~~~~~~~~~~CreateTree within LJ_pvfmm called~~~~~~~~~~~~~" << endl;

  //initialize the matricies, only needs to be done once
  pvfmm::PtFMM matrices(&mem_mgr);
  matrices.Initialize(mult_order, *comm, &LJ_ker);
  cout << "~~~~~~~~~~~~~~matrices.Initialize within LJ_pvfmm called~~~~~~~~~~~~~" << endl;

  //actually build tree
  tree->SetupFMM(&matrices);

  //evaluate tree
  PtFMM_Evaluate(tree, poten_and_grad, *N);
  cout << "~~~~~~~~~~~~~~LJ_pvfmm tree evaluated~~~~~~~~~~~~~" << endl;
  cout << "Values of first 4 numbers in poten_and_grad are: " << poten_and_grad[0] << ", " << poten_and_grad[1] << ", " << poten_and_grad[2] << " & " << poten_and_grad[3] << endl;

  //delete tree once finished
  delete tree;

  //##########

  //Comment all this out if you don't want to see this
  ofstream ffile;
  ffile.open("forces.txt");
  ffile << "fx, fy, fz (kg nm / ns^2)\n";
  for (i=0; i < *N ; i++)
    {
      (*force)[3*i  ] += poten_and_grad[4*i  ]         ;
      (*force)[3*i+1] += poten_and_grad[4*i+1]         ; 
      (*force)[3*i+2] += poten_and_grad[4*i+2]         ;
      *LJpoten_energy += poten_and_grad[4*i+3] * (*sig);                                                   //Final value of potential energy has a correction factor due to kernel definition
      ffile << (*force)[3*i]* *sl << "  " << (*force)[3*i+1]* *sl << "  " << (*force)[3*i+2]* *sl << "\n"; //All forces are re-scaled back to units of (N) by multiplying by side length 
    };

  return 0;
};
