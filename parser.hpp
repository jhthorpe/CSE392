// header file for our Parser class
#ifndef _PARSER_HPP_
#define _PARSER_HPP_

class Parser
{
  public:
    //num part, box length, temp, mass, timestep, num steps, LJ sigma, LJ epsilon, options
    int getInput(int *, float * ,float *, float *, float *, int*, float *, float *,int * );	//actual parser
};
#endif
