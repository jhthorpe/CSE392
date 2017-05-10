#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <mpi.h>
#include "verlet.hpp"
#include "potentials.hpp"
#include "energy.hpp"
#include "pvfmm.hpp"
using namespace std;

//Integrates position and velocity based on velocity verlet scheme
//For now, just for Lennard Jones and electrostatic forces between gas of atoms
//Requires prior initialization of 1D vectors mass, position and velocity 

void printer(){
  
};

void verlet:: Integration(int *N, int *nsteps, int *showsteps, double *tstep, double *sl, double *sig, double *eps, vector<double> *q, vector<double> *m, vector<double> *pos, vector<double> *vel, MPI_Comm* comm, int *start, int *end) {

  // Variables
  // N          : number of molecules local to this task 
  // nsteps     : number of time steps for integration
  // showsteps  : number of time steps for output display
  // tstep      : time step for integration
  // sl         : box length
  // sig        : sigma in Lennard Jones potential
  // eps        : epsilon in Lennard Jones potential
  // q          : charge on each atom (N)
  // m          : 1D vector of mass (size N)
  // pos        : 1D vector of positions (size 3N)
  // vel        : 1D vector of velocities (size 3N)
  // comm	: MPI communicator object
  // start	: starting index for this task
  // end	: ending index for this task

  //Get rank and maxrank
  int rank, size, maxrank;
  MPI_Comm_rank(*comm, &rank);
  MPI_Allreduce(&rank,&maxrank,1, MPI_INT, MPI_MAX, *comm);
  
  if (!rank) cout << "-------------Starting Verlet Integration-----------------" << endl;

  //Initial print a tops of files
  if (!rank){
      ofstream posfile, velfile, efile;
      posfile.open("traj_pos.xyz",ofstream::out);
      velfile.open("traj_vel.txt",ofstream::out);
      efile.open("traj_energy.txt",ofstream::out);
      posfile << "number of particles\n\n"; 
      velfile << "Time(ns) Atom xvel(nm/ns) yvel(nm/ns) zvel(nm/ns)\n";
      efile << "Time(ns) LJpot(kg-nm^2/ns^2) elcpot(kg-nm^2/ns^2) kin(kg-nm^2/ns^2) T(K)\n";
      posfile.close();
      velfile.close();
      efile.close();
  }

  //Start counting steps
  int step = 0;

  //Internal variables
  vector<double> f1(*N * 3), f2(*N * 3);    //Initialize temporary force arrays for each atomic coordinate   
  potentials pot;                           //potentials class object, pot 
  energy e;                                 //energy class object, e
  double kin, LJpot, elcpot, T;             //Declare energy variables

do {
  fill(f1.begin(),f1.end(),0.0); fill(f2.begin(),f2.end(),0.0);  //Initialize temporary force variables f1 and f2
  kin=0.0, LJpot=0.0, elcpot=0.0, T=0.0;                         //Re-initialize energy variables everytime

  pot.LJ_pvfmm(N,sl,sig,eps,pos,&f1,&LJpot,comm);     //Calculate LJ forces and potential energy on all N atoms at t=tstep*step
  pot.elc_pvfmm(N,sl,q,pos,&f1,&elcpot,comm,start);   //Calculate elc forces and potnetial energy on all N atoms at t=tstep*step
  e.kinetic_seq(N,sl,m,vel,&kin,&T);                  //Calculate total kinetic energy and temperature of system at t=tstep*step

  //Loop over all atoms to update position
  for (int ctr=0; ctr < *N * 3; ctr+=3) {

    //Update position of every atomic coordinate
    (*pos)[ctr] = (*pos)[ctr] + (*vel)[ctr]*(*tstep) + pow(*tstep,2)*f1[ctr]/(2 * (*m)[ctr/3]);
    (*pos)[ctr+1] = (*pos)[ctr+1] + (*vel)[ctr+1]*(*tstep) + pow(*tstep,2)*f1[ctr+1]/(2 * (*m)[ctr/3]);
    (*pos)[ctr+2] = (*pos)[ctr+2] + (*vel)[ctr+2]*(*tstep) + pow(*tstep,2)*f1[ctr+2]/(2 * (*m)[ctr/3]);

    //Correct for periodic boundaries
    //    (*pos)[ctr]-= floor((*pos)[ctr]/(*sl))*(*sl);
    //    (*pos)[ctr+1]-= floor((*pos)[ctr+1]/(*sl))*(*sl);
    //    (*pos)[ctr+2]-= floor((*pos)[ctr+2]/(*sl))*(*sl);
    (*pos)[ctr  ] -= floor((*pos)[ctr  ]);
    (*pos)[ctr+1] -= floor((*pos)[ctr+1]);
    (*pos)[ctr+2] -= floor((*pos)[ctr+2]);
  }

  pot.LJ_pvfmm(N,sl,sig,eps,pos,&f2,&LJpot,comm);      //Calculate LJ force on all N atoms after dt   
  pot.elc_pvfmm(N,sl,q,pos,&f2,&elcpot,comm,start);         //Calculate elc force on all N atoms after dt

  //Loop over all atoms to update velocity 
  for (int ctr=0; ctr < *N * 3; ctr+=3) {

    //Update velocity of every atomic coordinate
    (*vel)[ctr] = (*vel)[ctr] + (*tstep)*(f1[ctr]+f2[ctr])/(2 * (*m)[ctr/3]);
    (*vel)[ctr+1] = (*vel)[ctr+1] + (*tstep)*(f1[ctr+1]+f2[ctr+1])/(2 * (*m)[ctr/3]);
    (*vel)[ctr+2] = (*vel)[ctr+2] + (*tstep)*(f1[ctr+2]+f2[ctr+2])/(2 * (*m)[ctr/3]);

  }

  //check if need to print output
  if (step%(*showsteps)==0) {
    MPI_Barrier(*comm);		//this hurts my soul
    for(int j=0;j<=maxrank;j++){
      if (rank == j) {
        ofstream posfile, velfile, efile; 
        posfile.open("traj_pos.xyz", ofstream::out|ofstream::app);
        velfile.open("traj_vel.txt", ofstream::out|ofstream::app);
        efile.open("traj_energy.txt", ofstream::out|ofstream::app);
      
        for (int k=0; k < *N * 3; k+=3){
          posfile << "Ar " << (*pos)[k]*(*sl) << " " << (*pos)[k+1]*(*sl) << " " << (*pos)[k+2]*(*sl) << "\n";
          velfile << (*tstep)*step << " " << k/3 << " " << (*vel)[k]*(*sl) << " " << (*vel)[k+1]*(*sl) << " " << (*vel)[k+2]*(*sl) << "\n";
          efile << (*tstep)*step << " " << LJpot*(*sl * *sl) << " " << elcpot*(*sl * *sl) << " " << kin*(*sl * *sl) << " " << T << "\n";  //Rescale output to correct units (Joules)
        }
        posfile.close();
        velfile.close();
        efile.close();
      }
      MPI_Barrier(*comm);
    } 
  }
  // I am fully aware that this risks killing a fair amount of our parallel scaling.
   

  //Update time step
  step+=1;


 } while(step < *nsteps);//End do-while loop after nsteps

if (!rank) cout << "----------Dynamics Ended----------------" << endl;


}
