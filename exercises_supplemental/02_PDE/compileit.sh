#!/bin/sh

# optimized build
#g++ -O3 -DNDEBUG -o predator-prey-pde predator-prey-pde.cpp -lboost_program_options

# debug build
g++ -O0 -g -o predator-prey-pde predator-prey-pde.cpp -lboost_program_options