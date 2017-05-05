//header file for forces class
#ifndef _FORCES_HPP_
#define _FORCES_HPP_
#include <vector>
using namespace std;

class Forces
{
  public:
    //sequential LJ forces
    //number of mol, box length, sigma param, epsilon param, position vector, forces vector 
    int LJ_seq(int *, double *, double *, double *, vector<double> * , vector<double> *); 

    //sequential coulombic intecactions
    //number of mol, box length, charges, position vector, forces vector 
    int elc_seq(int *,double *, vector<double> *, vector<double> *, vector<double> *);

    //sequential LJ forces with periodic boundaries
    //number of mol, box length, sigma param, epsilon param, position vector, forces vector 
    int LJ_seq_bound(int *, double *, double *, double *, vector<double> * , vector<double> *);

    //sequential coulombic intecactions
    //number of mol, box length, charges, position vector, forces vector 
    int elc_seq_bound(int *,double *, vector<double> *, vector<double> *, vector<double> *);

    //parallel LJ forces with periodic boundaries
    //number of mol, box length, sigma param, epsilon param, position vector, forces vector 
    int LJ_omp_bound(int *, double *, double *, double *, vector<double> * , vector<double> *);
};


#endif
