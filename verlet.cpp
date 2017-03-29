#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;

//Define function force for particle
double force(double x) {
  return -x;
}

//Define function Verlet to compute dynamics
void Verlet(vector<double> mass, vector<double> position, vector<double> velocity, double tstep, double nsteps) {

  //Open files for storing position and velocity
  ofstream posfile, velfile;
  posfile.open("traj_pos.txt"); velfile.open("traj_vel.txt");

  //Start counting steps
  int step = 0;
  cout << "Following are the coordinates after every 10 time steps: " << endl;
  //Velocity Verlet
  do {

    //Loop over all atoms
    for (int ctr=0; ctr<position.size(); ctr++) {

      //Update position
      double f1 = force(position[ctr]);        
      position[ctr] = position[ctr] + velocity[ctr]*tstep + pow(tstep,2)*f1/(2*mass[ctr/3]);

      //Update velocity
      double f2 = force(position[ctr]);
      velocity[ctr] = velocity[ctr] + tstep*(f1+f2)/(2*mass[ctr/3]);

      //Write every ten steps to file
      if (step%10==0) {
	posfile<< step*tstep << " " << position[ctr] << "\n";
	velfile<< step*tstep << " " << velocity[ctr] << "\n";
	cout << "t=" << step*tstep << ": x[" << ctr << "]=" << position[ctr] << ", v[" << ctr << "]=" << velocity[ctr] << endl;
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
  Verlet(mass, x, v, tstep, nsteps);

  return 0;
}
