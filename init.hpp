// header file for our Init (initializer) class
#ifndef _INIT_HPP_
#define _INIT_HPP_
#include <vector>
using namespace std;

class Init
{
  public:
    // number of mol, box length, temp, position vector, velocity vector
    int initialize(int *, float *, float *, vector<float> *, vector<float> * );
};

#endif

