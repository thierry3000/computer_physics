#include<stdlib.h>
#include<string.h>
#include<stdio.h>
#include<limits.h>
#include<time.h>
#include<cmath>
#include<assert.h>
#include<iostream>
#include<map>
#include<fstream>
#include<iomanip>
#include"nr/nr.h"
#include"nr/ran3.cpp"
#include <cstdlib>
#include <cstdio>




using namespace std;

int main() 
{
 char dateiname[50]="Ausgabe.dat";
 fstream Datei;
 Datei.open(dateiname, ios::out);

 unsigned long long i, nr_samples;

 nr_samples=100;
 int seed=1;
 
 for(i=0;i<nr_samples;i++)
 {
  Datei<<NR::ran3(seed)<<endl;

 }
 Datei.close();

}


