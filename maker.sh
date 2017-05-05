#! /bin/bash
mpicxx -o runMD -std=c++11 runMD.cpp parser.cpp killer.cpp init.cpp forces.cpp -fopenmp energy.cpp verlet.cpp forces_pvfmm.cpp 
