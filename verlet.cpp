#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
//Developer defined header files
#include "verlet.hpp"
#include "forces.hpp"
#include "energy.hpp"
using namespace std;

//Integrates position and velocity based on velocity verlet scheme
//For now, just for Lennard Jones and electrostatic forces between gas of atoms
//Requires prior initialization of 1D vectors mass, position and velocity 

void verlet:: Integration(int *N, int *nsteps, double *tstep, double *sl, double *sig, double *eps, vector<double> *q, vector<double> *m, vector<double> *pos, vector<double> *vel) {

  // Varaibles                                                                                                                                          
  // N          : number of molecules                                                                                                                                        
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

  //Open files for storing position and velocity
  ofstream posfile, velfile, kefile, pefile;
  posfile.open("traj_pos.txt"); velfile.open("traj_vel.txt"); kefile.open("traj_kenergy.txt"); pefile.open("traj_penergy.txt");

  posfile << "Time(ns) Atom xpos(nm) ypos(nm) zpos(nm)\n"; velfile << "Time(ns) Atom xvel(nm/ns) yvel(nm/ns) zvel(nm/ns)\n"; 
  kefile << "Time(ns) T(K) kin(kg-nm^2/ns^2)\n"; pefile << "Time(ns) LJpot(kg-nm^2/ns^2) elcpot(kg-nm^2/ns^2)\n";

  cout << "-------------Starting Verlet Integration-----------------" << endl;

  //Start counting steps
  int step = 0;

  //Internal variables
  int showsteps=1000;
  vector<double> f1(*N * 3), f2(*N * 3);   //Initialize temporary force arrays for each atomic coordinate   
  Forces forces;                           //Forces class object, forces 
  energy e;                                //energy class object, e
  double kin, LJpot, elcpot, T;            //Declare energy variables

do {
  fill(f1.begin(),f1.end(),0.0); fill(f2.begin(),f2.end(),0.0);  //Initialize temporary force variables f1 and f2
  kin=0.0, LJpot=0.0, elcpot=0.0, T=0.0;                         //Re-initialize energy variables everytime

  forces.LJ_seq_bound(N,sl,sig,eps,pos,&f1);      //Calculate LJ force on all N atoms at t=tstep*step
  forces.elc_seq_bound(N,sl,q,pos,&f1);           //Calculate elc force on all N atoms at t=tstep*step

  //Write energies to file every showsteps steps
  if (step%showsteps==0) {

    e.LJpot_seq(N,sl,sig,eps,pos,&LJpot);     //Calculate total lennard jones potential energy of system at t=tstep*step                                                           
    e.elcpot_seq(N,q,pos,&elcpot);            //Calculate total electrostatic potential energy of system at t=tstep*step                                                           
    e.kinetic_seq(N,m,vel,&kin,&T);           //Calculate total kinetic energy and temperature of system at t=tstep*step 

    pefile << (*tstep)*step << " " << LJpot << " " << elcpot << "\n";  
    kefile << (*tstep)*step << " " << T << " " << kin << "\n";
  }

  //Loop over all atoms to update position
  for (int ctr=0; ctr < *N * 3; ctr+=3) {

    //Write every showsteps steps to file
    if (step%showsteps==0) {
      posfile << (*tstep)*step << " " << ctr/3 << " " << (*pos)[ctr] << " " << (*pos)[ctr+1] << " " << (*pos)[ctr+2] << "\n";
      cout << "t = " << (*tstep)*step*1000 << " ps" << ", Atom number " << ctr/3 << ": x=" << (*pos)[ctr] << ", y=" << (*pos)[ctr+1] << ", z=" << (*pos)[ctr+2] << endl;
    }//End if statement for writing to file                                                                                                                                        

    //Update position of every atomic coordinate
    (*pos)[ctr] = (*pos)[ctr] + (*vel)[ctr]*(*tstep) + pow(*tstep,2)*f1[ctr]/(2 * (*m)[ctr/3]);
    (*pos)[ctr+1] = (*pos)[ctr+1] + (*vel)[ctr+1]*(*tstep) + pow(*tstep,2)*f1[ctr+1]/(2 * (*m)[ctr/3]);
    (*pos)[ctr+2] = (*pos)[ctr+2] + (*vel)[ctr+2]*(*tstep) + pow(*tstep,2)*f1[ctr+2]/(2 * (*m)[ctr/3]);

    //Correct for periodic boundaries                                                                                                                                              
    (*pos)[ctr]-= floor((*pos)[ctr]/(*sl))*(*sl);
    (*pos)[ctr+1]-= floor((*pos)[ctr+1]/(*sl))*(*sl);
    (*pos)[ctr+2]-= floor((*pos)[ctr+2]/(*sl))*(*sl);

  }//End for loop over all particles

  forces.LJ_seq_bound(N,sl,sig,eps,pos,&f2);      //Calculate LJ force on all N atoms after dt   
  forces.elc_seq_bound(N,sl,q,pos,&f2);           //Calculate elc force on all N atoms after dt

  //Loop over all atoms to update velocity 
  for (int ctr=0; ctr < *N * 3; ctr+=3) {

    //Write every showsteps steps to file
    if (step%showsteps==0) {
      velfile << (*tstep)*step << " " << ctr/3 << " " << (*vel)[ctr] << " " << (*vel)[ctr+1] << " " << (*vel)[ctr+2] << "\n";
      cout << "t = " << (*tstep)*step*1000 << " ps" << ", Atom number " << ctr/3 << ": vx=" << (*vel)[ctr] << ", vy=" << (*vel)[ctr+1] << ", vz=" << (*vel)[ctr+2] << endl;
    }//End if statement for writing to file 

    //Update velocity of every atomic coordinate
    (*vel)[ctr] = (*vel)[ctr] + (*tstep)*(f1[ctr]+f2[ctr])/(2 * (*m)[ctr/3]);
    (*vel)[ctr+1] = (*vel)[ctr+1] + (*tstep)*(f1[ctr+1]+f2[ctr+1])/(2 * (*m)[ctr/3]);
    (*vel)[ctr+2] = (*vel)[ctr+2] + (*tstep)*(f1[ctr+2]+f2[ctr+2])/(2 * (*m)[ctr/3]);

  }//End for loop over all particles                                                                                                                                               

  //Update time step                                                                                                                                                               
  step+=1;

 } while(step < *nsteps);//End do-while loop after nsteps                                                                                                                          

 cout << "----------Dynamics Ended----------------" << endl;

//Close file streams                                                                                                                                                            
 posfile.close(); velfile.close(); kefile.close(); pefile.close();

}
