//header file for forces class
#ifndef _FORCES_HPP_
#define _FORCES_HPP_
#include <vector>
using namespace std;

class Forces
{
  public:
    //input stuff
    void Test();

    //sequential LJ forces
    //number of mol, box length, sigma param, epsilon param, position vector, forces vector 
    int LJ_seq(int *, double *, double *, double *, vector<double> * , vector<double> *); 
};

#endif
