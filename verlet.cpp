#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;

int boxx=10, boxy=10, boxz=10;   //Box Size
double cutoff=3.0, eps=1.0, sgm=1.0;    // LennardJones
vector<double> mass(2,1.0),position={3.0,5.0,5.0,7.0,5.0,5.0},velocity(6,0.0);      //Mass, Positions, Velocities

//LennardJones interaction between two atoms
vector<double> vdw(double x, double y, double z, double x1, double y1, double z1) {
  vector<double> vander(3,0.0);
  double rx=x-x1,ry=y-y1,rz=z-z1;
  double r=sqrt(rx*rx+ry*ry+rz*rz);
  if (r<cutoff) {
    double f=-(24*eps/r)*(2*pow(sgm/r,12)-pow(sgm/r,6));
    vander={f*(rx/r),f*(ry/r),f*(rz/r)};
  }
  return vander;
}

//Mean Force on particle
vector<double> force(double x, double y, double z) {
  vector<double> f(3,0.0), LJ(3);
  for (int ctr=0; ctr<position.size();ctr=ctr+3) {
    if ((x==position[ctr])&&(y==position[ctr+1])&&(z==position[ctr+2])) {
      continue;
    }
    else {
      LJ=vdw(x,y,z,position[ctr],position[ctr+1],position[ctr+2]);
      for (int i=0;i<3;i++) {
        f[i]+=LJ[i];
      }
    }
  }
  return f;
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

      position[ctr]-= floor(position[ctr]/boxx)*boxx; 
      position[ctr+1]-= floor(position[ctr+1]/boxy)*boxy;
      position[ctr+2]-= floor(position[ctr+2]/boxz)*boxz;
       
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

  double tstep=0.01;                                 //Define time step                                                                                       
  double nsteps=100;                                 //Define number of steps
  int key;
  cout << "Press Enter to start simulation: ";
  cin >> key;
  //Run Verlet function
  Verlet(mass, position, velocity, tstep, nsteps, boxx, boxy, boxz);

  return 0;
}
