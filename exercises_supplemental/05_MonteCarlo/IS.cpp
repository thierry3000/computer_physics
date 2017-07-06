/*
 * part of computer physics lecture at Saarland University
 *                                  by Prof. Dr. Heiko Rieger and Adam Wysocki
 * 
 *      summer term 2017
 * 
 * exercises by Thierry Fredrich
 * 
 * 
 * compile with:
 * g++ IS.cpp -std=c++0x -o myIsingSimulation -lboost_program_options
 */

#include<stdlib.h>
#include<string.h>
#include<stdio.h>
#include<assert.h>
#include<iostream>
#include <vector>

// since the compiler in the cip pool is rather old
// we need some hacking.
#include<features.h> // gives various information on the build system
#if __GNUC_PREREQ(4,6)
    // means gnu compiler is higher than 4.5 and fully c++11 compatible
    // we don't need to hacking
#else
  //old fashioned
  #define OLD_RANDOM
  #include<cstdlib>
  #include<time.h>
#endif
// c++ version von printf.
#include <boost/format.hpp>
// Boost program_options support command line arguments.
#include <boost/program_options.hpp>

//random numbers
/* there are tons of random number generators!
 * in fact this is a research topic on its own.
 * Here: we keep it simple for now.
 */
#include <random>
/* note:
 * it is not the fastest way to write this in a external function,
 * but it keeps the focus on the physical stuff.
 */
double create_random_number()
{
#ifdef OLD_RANDOM
  std::srand(std::time(0));
  return (double)(std::rand()/RAND_MAX);
#else
  //from http://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0, 1);
  return dis(gen);
#endif
}

//declaration
typedef std::vector<bool> zustand;

using namespace std;

int seed=3;
int a=8;
int b=8;
double J=1;
double beta=1;
double p=1;
int warmlauf=10000;

// this is a version where no boost library is used for pasing the command line arguments. Not that this could cause serious trouble since arguments are not checked!
void ReadCommandLine(int argc, char** argv, int & steps){
  for( int i = 1; i<argc; i++ ){
    if (strstr(argv[i], "-la" )) {a=atoi(argv[i]+3);}
    if (strstr(argv[i], "-lb" )) {b=atoi(argv[i]+3);}
    if (strstr(argv[i], "-seed" )) {seed=atoi(argv[i]+5);}
    if (strstr(argv[i], "-steps" )) {steps=atoi(argv[i]+6);}
    if (strstr(argv[i], "-beta" )) {beta=atof(argv[i]+5);}
    if (strstr(argv[i], "-p" )) {p=atof(argv[i]+2);}
    if (strstr(argv[i], "-warm" )) {warmlauf=atoi(argv[i]+5);}
  }
}

void Initialisierung(zustand &z, double p){
  for(int i=0; i<a*b; i++){
    if((create_random_number())<p){
      z[i]=true;
    }else{
      z[i]=false;
    }
  }
}

double dE(const zustand & z, int i){
  int sum;
  sum=(z[i]==z[(i+1)-(((i+1)%b)==0)*b])+(z[i]==z[(i-1)+((i%b)==0)*b])+(z[i]==z[(i+a)%(a*b)])+(z[i]==z[(i-a+a*b)%(a*b)]);
  return -2*(4.0-2.0*sum);
}

zustand Metropolis_step(zustand & z){
  int turn=floor((create_random_number())*a*b);
  if (turn==a*b){
    --turn;
  }
  double za=(create_random_number());
  if (exp(-beta*dE(z, turn))>=za){
    z[turn]=!z[turn];
  }
  return z;
}

int Magnetisierung(const zustand & z){  
  int counter=0;
  for(int i=0;i<a*b;++i){
    if (z[i]){
      ++counter;
    }else{
      --counter;
    }
  }
  return counter;
}

/* exercise 2. a) */
// double Magnetisierung(const zustand & z){
//   return 42.0;
// }

int main(int argc, char** argv){
  int steps=1;
  double dmaggi=0.0;
  ReadCommandLine(argc, argv, steps);
  vector<bool> z(a*b);
  vector<int> m(steps+1);
  Initialisierung(z, p);
  m[0]=Magnetisierung(z);

  for(int counter=1; counter < warmlauf+1; ++counter){
    z=Metropolis_step(z);
    if(counter%100==0){
      cout << "warmlauf:" << counter << endl;
    }
  }

  for(int counter=1; counter<steps+1; ++counter){
    z=Metropolis_step(z);
    dmaggi+=abs(m[counter]=Magnetisierung(z));
  }

  cout << dmaggi/(steps*a*b) << endl;
}

