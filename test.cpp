#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

double cutoff=10.0, eps=1.0, sgm=1.0;

vector<double> vdw(double x, double y, double z, double x1, double y1, double z1) {   //LennardJones interaction between two atoms                                                   
  vector<double> vander(3,0.0); double f;
  //Find relative distances                                                                      
  double rx=x-x1,ry=y-y1,rz=z-z1;
  double r=sqrt(rx*rx+ry*ry+rz*rz);
  //Truncate effect of Lennard-Jones after certain cut-off                                       
  if (r<cutoff) {
    f = (24*eps/r)*(2*pow(sgm/r,12)-pow(sgm/r,6));
    vander={f*(rx/r),f*(ry/r),f*(rz/r)};
  }//End if statement                                                                                                                                                                

  return vander;
}//End function vdw()

int main() {
  vector<double > x,x1,LJ;
  x={3.0,5.0,5.0};x1={7.0,5.0,5.0};
  cout << "Vector x has following elements: " << endl;
  for (int i=0; i<x.size(); i++) {
    cout << x[i] << endl;
  }
  //double x=-0.2, y=1.0;
  //cout << "x/y " << x/y << endl;
  //cout << "(floor)x/y " << floor(x/y) << endl;
  cout << "Now calling Verlet()... " << endl;
  LJ = vdw(x[0],x[1],x[2],x1[0],x1[1],x1[2]);
  cout << "Force vector LJ has following elements:" << endl;
  for (int j=0; j<LJ.size(); j++) {
    cout << LJ[j] << endl;
  }
  return 0;
}
