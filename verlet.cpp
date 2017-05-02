#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
//Developer defined header files
#include "verlet.hpp"
#include "forces.hpp"
using namespace std;

//Integrates position and velocity based on velocity verlet scheme
//For now, just for Lennard Jones force between gas of neutral atoms
//Requires prior initialization of 1D vectors mass, position and velocity 

void verlet:: Integration(int *N, int *nsteps, double *tstep, double *boxL, double *sig, double *eps, vector<double> *m, vector<double> *pos, vector<double> *vel) {

  // Varaibles                                                                                                                                          
  // N          : number of molecules                                                                                                                                        
  // nsteps     : number of time steps for integration
  // tstep      : time step for integration
  // boxL       : box length
  // sig        : sigma in Lennard Jones potential
  // eps        : epsilon in Lennard Jones potential
  // m          : 1D vector of mass (size N)                                                                                                                                 
  // pos        : 1D vector of positions (size 3N)                                                                                                                    
  // vel        : 1D vector of velocities (size 3N)                                                                                            

  //Open files for storing position and velocity
  ofstream posfile, velfile;
  posfile.open("traj_pos.txt"); velfile.open("traj_vel.txt");

  posfile << "Time xpos ypos zpos\n"; velfile << "Time xvel yvel zvel\n"; 

  cout << "-------------Starting Verlet Integration-----------------" << endl;

  //Start counting steps
  int step = 0;
                                                                                                                                                                  
do {

  //Initialize temporary force arrays for each atomic coordinate
  vector<double> f1(*N * 3, 0.0), f2(*N * 3, 0.0);

  Forces forces;        //Forces class object, forces

  forces.LJ_seq_bound(N,boxL,sig,eps,pos,&f1);      //Calculate force on all N atoms at t=tstep*step

  //Loop over all atoms to update position
  for (int ctr=0; ctr < *N * 3; ctr+=3) {

    //Write every ten steps to file
    if (step%10==0) {
      posfile << (*tstep)*step << " " << (*pos)[ctr] << " " << (*pos)[ctr+1] << " " << (*pos)[ctr+2] << "\n";
      cout << " Following are the coordinates after " << (*tstep)*step << "femtoseconds" << endl;
      cout << " x=" << (*pos)[ctr] << ", y=" << (*pos)[ctr+1] << ", z=" << (*pos)[ctr+2] << endl;
    }//End if statement for writing to file                                                                                                                                        

    //Update position of every atomic coordinate
    (*pos)[ctr] = (*pos)[ctr] + (*vel)[ctr]*(*tstep) + pow(*tstep,2)*f1[ctr]/(2 * (*m)[ctr/3]);
    (*pos)[ctr+1] = (*pos)[ctr+1] + (*vel)[ctr+1]*(*tstep) + pow(*tstep,2)*f1[ctr+1]/(2 * (*m)[ctr/3]);
    (*pos)[ctr+2] = (*pos)[ctr+2] + (*vel)[ctr+2]*(*tstep) + pow(*tstep,2)*f1[ctr+2]/(2 * (*m)[ctr/3]);

    //Correct for periodic boundaries                                                                                                                                              
    (*pos)[ctr]-= floor((*pos)[ctr]/(*boxL))*(*boxL);
    (*pos)[ctr+1]-= floor((*pos)[ctr+1]/(*boxL))*(*boxL);
    (*pos)[ctr+2]-= floor((*pos)[ctr+2]/(*boxL))*(*boxL);

  }//End for loop over all particles

  forces.LJ_seq_bound(N,boxL,sig,eps,pos,&f2);      //Calculate force on all N atoms after dt   

  cout << "-------------------------------------------------------------" << endl;

  //Loop over all atoms to update velocity 
  for (int ctr=0; ctr < *N * 3; ctr+=3) {

    //Write every ten steps to file
    if (step%10==0) {
      velfile << (*tstep)*step << " " << (*vel)[ctr] << " " << (*vel)[ctr+1] << " " << (*vel)[ctr+2] << "\n";
      cout << " Following are the velocities after " << (*tstep)*step << "femtoseconds" << endl;
      cout << " vx=" << (*vel)[ctr] << ", vy=" << (*vel)[ctr+1] << ", vz=" << (*vel)[ctr+2] << endl;
    }//End if statement for writing to file 

    //Update velocity of every atomic coordinate
    (*vel)[ctr] = (*vel)[ctr] + (*tstep)*(f1[ctr]+f2[ctr])/(2 * (*m)[ctr/3]);
    (*vel)[ctr+1] = (*vel)[ctr+1] + (*tstep)*(f1[ctr+1]+f2[ctr+1])/(2 * (*m)[ctr/3]);
    (*vel)[ctr+2] = (*vel)[ctr+2] + (*tstep)*(f1[ctr+2]+f2[ctr+2])/(2 * (*m)[ctr/3]);

  }//End for loop over all particles                                                                                                                                               

  //Update time step                                                                                                                                                               
  step+=1;

 } while(step < *nsteps);//End do-while loop after nsteps                                                                                                                          

//Close file streams                                                                                                                                                            
 posfile.close(); velfile.close();

}
