#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;

double boxx=10.0, boxy=10.0, boxz=10.0;   //Box Size (nm)

double cutoff=10.0, eps=0.01, sgm=0.3;    // LennardJones (Arbitrary units for now)

vector<double> mass={1.0,10.0}, position={1.5,2.4,3.8,4.9,5.2,6.1}, velocity(6,0.0);      //Mass, Positions, Velocities
                                                                                      
vector<double> vdw(double x, double y, double z, double x1, double y1, double z1) {   //LennardJones interaction between two atoms
  vector<double> vander(3,0.0); double f;
                                                                                      //Find relative distances
  double rx=x-x1,ry=y-y1,rz=z-z1, r, diag=sqrt(boxx*boxx+boxy*boxy+boxz*boxz);

  //Correct for boundary conditions
  do {
  rx-= floor(rx/boxx)*boxx;
  ry-= floor(ry/boxy)*boxy;
  rz-= floor(rz/boxz)*boxz;
  r=sqrt(rx*rx+ry*ry+rz*rz);
  } while(r>diag); 
                                                                                     //Truncate effect of Lennard-Jones after certain cut-off
  if (r<cutoff) {
    f = (24*eps/r)*(2*pow(sgm/r,12)-pow(sgm/r,6));
    vander={f*(rx/r),f*(ry/r),f*(rz/r)};
  }//End if statement

  return vander;
}//End function vdw()

vector<double> force(double x, double y, double z) {                                  //Mean Force on particle
  vector<double> f(3,0.0), LJ(3);

  for (int ctr=0; ctr<position.size();ctr=ctr+3) {                                    //Calculate forces from all particles

    if (x!=position[ctr]) {                                                           //Skip calculation of force due to self

      LJ=vdw(x,y,z,position[ctr],position[ctr+1],position[ctr+2]);                    //Store LJ force due to atom number 'ctr/3' 

      for (int i=0;i<3;i++) {                                                         //Add up all forces on each component 
        f[i]+=LJ[i];
      }//End for loop over force components
    }//End if statement
  }//End for loop over all particles
  return f;
}//End function force()

//Define function Verlet to compute dynamics
void Verlet(vector<double> mass, vector<double> position, vector<double> velocity, double tstep, double nsteps, double boxx, double boxy, double boxz) {

  cout << "Verlet started successfully..." << endl;

  //Temperory variables for forces
  vector<double> f1(3),f2(3);

  //Open files for storing position and velocity
  ofstream posfile, velfile;
  posfile.open("traj_pos.txt"); velfile.open("traj_vel.txt");

  //Start counting steps
  int step = 0;
  cout << "Following are the coordinates after every time step: " << endl;

  //Velocity Verlet
  do {

    //Loop over all atoms
    for (int ctr=0; ctr<position.size(); ctr=ctr+3) {

      //Write every ten steps to file                                                                                                                                                
      if (step%10==0) {
        posfile<< step*tstep << " " << position[ctr] << position[ctr+1] << position[ctr+2] << "\n";
        velfile<< step*tstep << " " << velocity[ctr] << velocity[ctr+1] << velocity[ctr+2] << "\n";
        cout << "t=" << step*tstep;
        cout << ": x=" << position[ctr] << ", y=" << position[ctr+1] << ", z=" << position[ctr+2] << endl;
        cout << " vx=" << velocity[ctr] << ", vy=" << velocity[ctr+1] << ", vz=" << velocity[ctr+2] << endl;
      }//End if statement for writing to file

      //Update position
      f1 = force(position[ctr], position[ctr+1], position[ctr+2]);        
      position[ctr] = position[ctr] + velocity[ctr]*tstep + pow(tstep,2)*f1[0]/(2*mass[ctr/3]);
      position[ctr+1] = position[ctr+1] + velocity[ctr+1]*tstep + pow(tstep,2)*f1[1]/(2*mass[ctr/3]);
      position[ctr+2] = position[ctr+2] + velocity[ctr+2]*tstep + pow(tstep,2)*f1[2]/(2*mass[ctr/3]);

      //Correct for periodic boundaries
      position[ctr]-= floor(position[ctr]/boxx)*boxx; 
      position[ctr+1]-= floor(position[ctr+1]/boxy)*boxy;
      position[ctr+2]-= floor(position[ctr+2]/boxz)*boxz;
       
      //Update velocity
      f2 = force(position[ctr], position[ctr+1], position[ctr+2]);
      velocity[ctr] = velocity[ctr] + tstep*(f1[0]+f2[0])/(2*mass[ctr/3]);
      velocity[ctr+1] = velocity[ctr+1] + tstep*(f1[1]+f2[1])/(2*mass[ctr/3]);
      velocity[ctr+2] = velocity[ctr+2] + tstep*(f1[2]+f2[2])/(2*mass[ctr/3]);

    }//End for loop over all particles

    //Update time step
    step+=1;

  } while(step<nsteps);//End do-while loop over time

  //Close files
  posfile.close(); velfile.close();

}//End function Verlet()

//Dummy main function to test Verlet
int main() {

  //  vector<double> mass(1,1.0), position={5.0,5.0,5.0}, velocity(3,0.0);      //Mass, Positions, Velocities  
  double boxx=10.0, boxy=10.0, boxz=10.0;   //Box Size (nm)
  double tstep=0.000001;                        //Define time step (ns)
  int nsteps=100;                         //Define number of steps

  //  int key;
  //cout << "Press Enter to start simulation: ";
  //cin >> key;

  //Run Verlet function
  Verlet(mass, position, velocity, tstep, nsteps, boxx, boxy, boxz);

  return 0;
}//End function main()
