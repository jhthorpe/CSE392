#! /bin/bash
g++ -o runMD -std=c++11 runMD.cpp parser.cpp killer.cpp init.cpp forces.cpp -fopenmp energy.cpp verlet.cpp
