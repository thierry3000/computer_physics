cmake_minimum_required(VERSION 2.6)
project(comphys_supplemental)
set(CMAKE_CXX_STANDARD 11) # C++11...

add_subdirectory(02_PDE/)
add_subdirectory(04_Sampling/)
add_subdirectory(05_MonteCarlo/)
add_subdirectory(06_MolecularDynamics/)
add_subdirectory(07_QuantumMonteCarlo/)
add_subdirectory(08_ExactDiagonalization/)
