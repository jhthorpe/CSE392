#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;

int boxx=1.0, boxy=2.0, boxz=3.0;

//Define function force for particle
vector<double> electroforce(double x, double y, double z) {
  vector<double> electro;

  return electro;
}

//Define function Verlet to compute dynamics
void Verlet(vector<double> mass, vector<double> position, vector<double> velocity, double tstep, double nsteps, int boxx, int boxy, int boxz) {

  //Open files for storing position and velocity
  ofstream posfile, velfile;
  posfile.open("traj_pos.txt"); velfile.open("traj_vel.txt");

  //Start counting steps
  int step = 0;
  cout << "Following are the coordinates after every 10 time steps: " << endl;
  //Velocity Verlet
  do {

    //Loop over all atoms
    for (int ctr=0; ctr<position.size(); ctr=ctr+3) {

      //Update position
      double f1 = force(position[ctr], position[ctr+1], position[ctr+2]);        
      position[ctr] = position[ctr] + velocity[ctr]*tstep + pow(tstep,2)*f1[0]/(2*mass[ctr/3]);
      position[ctr+1] = position[ctr+1] + velocity[ctr+1]*tstep + pow(tstep,2)*f1[1]/(2*mass[ctr/3]);
      position[ctr+2] = position[ctr+2] + velocity[ctr+2]*tstep + pow(tstep,2)*f1[2]/(2*mass[ctr/3]);

      position[ctr]+= floor(position[ctr]/boxx)*boxx; 
      position[ctr+1]+= floor(position[ctr+1]/boxy)*boxy;
      position[ctr+2]+= floor(position[ctr+2]/boxz)*boxz;
       
      //Update velocity
      double f2 = force(position[ctr], position[ctr+1], position[ctr+2]);
      velocity[ctr] = velocity[ctr] + tstep*(f1[0]+f2[0])/(2*mass[ctr/3]);
      velocity[ctr+1] = velocity[ctr+1] + tstep*(f1[1]+f2[1])/(2*mass[ctr/3]);
      velocity[ctr+2] = velocity[ctr+2] + tstep*(f1[2]+f2[2])/(2*mass[ctr/3]);

      //Write every ten steps to file
      if (step%10==0) {
	posfile<< step*tstep << " " << position[ctr] << position[ctr+1] << position[ctr+2] << "\n";
	velfile<< step*tstep << " " << velocity[ctr] << velocity[ctr+1] << velocity[ctr+2] << "\n";
	cout << "t=" << step*tstep;
	cout << ": x=" << position[ctr] << ", y=" << position[ctr+1] << ", z=" << position[ctr+2] << endl;
	cout << " vx=" << velocity[ctr] << ", vy=" << velocity[ctr+1] << ", vz=" << velocity[ctr+2] << endl;
      }
    }

    //Update time step
    step+=1;
  } while(step<nsteps);

  //Close files
  posfile.close(); velfile.close();
}

//Dummy main function to test Verlet
int main() {

  vector<double> mass(1,1.0),x(3,0.0),v(3,1.0);      //Define mass, initial position and velocity                                                                                      
  double tstep=0.01;                                 //Define time step                                                                                       
  double nsteps=100;                                 //Define number of steps
  int key;
  cout << "Press Enter to start simulation: ";
  cin >> key;
  //Run Verlet function
  Verlet(mass, x, v, tstep, nsteps, boxx, boxy, boxz);

  return 0;
}
