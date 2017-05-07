// header file for our Init (initializer) class
#ifndef _INIT_HPP_
#define _INIT_HPP_
#include <vector>
using namespace std;

class Init
{
  public:
    // number of mol, box length, temp, mass, position vector, velocity vector
    int initialize(int *, double *,double *, vector<double> *, vector<double> *, vector<double> * );
};

#endif
