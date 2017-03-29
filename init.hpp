// header file for our Init (initializer) class
#ifndef _INIT_HPP_
#define _INIT_HPP_

class Init
{
  public:
    // number of mol, box length, temp, position vector, velocity vector
    int initialize(int *, float *, float *, float *, float *);
    int tester(float *);
};

#endif

