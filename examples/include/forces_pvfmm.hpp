//header file for forces class
#ifndef _FORCES_PVFMM_HPP_
#define _FORCES_PVFMM_HPP_
#include <vector>
using namespace std;

class Forces_pvfmm
{
  public:
    //using pvfmm to impliment electrostatics
    //number of mol, box length, charges, position vector, forces vector 
    int elc_pvfmm(int *,double *, vector<double> *, vector<double> *, vector<double> *);
};


#endif
