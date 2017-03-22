# CSE392
Implementing a basic MD simulation that utilizes the pvfmm program to evaluate coulomb interactions

Basic Library Structure

1. INPUT/setup
  1.0 Input parser : parser.cpp
  1.1 construct the simulation space (initial variables/arrays) : constructor.cpp
  1.11a if random setup selected: distribute, randomly molecules in this simulation space
  1.11b if already have simulation space constructed, read in values
  
2. DOSTUFF

3. OUTPUT
  3.0 Output the simulation space, if we want to use this again later
  3.1 Output the final energy of the system
  3.2 Tracker file, which records the state of the system every x timesteps
  3.3 Video file?
