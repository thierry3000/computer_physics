//Kompilieren g++ -O3 Lanczos.cpp -o Lanczos

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<fstream>
#include<sstream>
#include<iostream>
#include<iomanip>
#include<cmath>
#include<string>
#include<cstdlib>
#include<vector>
#include<boost/numeric/ublas/matrix.hpp>
#include<boost/numeric/ublas/io.hpp>
#include<boost/format.hpp>

// since the compiler in the cip pool is rather old
// we need some hacking.
#include<features.h> // gives various information on the build system
#if __GNUC_PREREQ(4,7)
//#if 0
  //... code requiring gcc 4.7 or later ...
    // means gnu compiler is higher than 4.5 and fully c++11 compatible
    // we don't need to hacking
  #include <random>
#else
  #define OLD_RANDOM
  // get random numbers from boost -- > cip pool is boost 1.41!!!!!
  #include <boost/random.hpp>
#endif

//random numbers
/* there are tons of random number generators!
 * in fact this is a research topic on its own.
 * Here: we keep it simple for now.
 */
/* note:
 * it is not the fastest way to write this in a external function,
 * but it keeps the focus on the physical stuff.
 */
double create_random_number()
{
#ifdef OLD_RANDOM  // we use boost to create random numbers
  /** see:
   * http://www.bnikolic.co.uk/blog/cpp-boost-uniform01.html
   */
  boost::mt19937 rng(std::time(0));
  static boost::uniform_01<boost::mt19937> uni(rng);
  return uni();
#else
  /** meanwhile this is include in the standard libraries
   * from http://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
   */
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0, 1);
  return dis(gen);
#endif
}

// note: for solving Tridiagonalmatrix see Thomas Algorithm
// here: some algrithms from numerical recepies, see book
namespace NR
{
  template<class T>
  inline const T SIGN(const T &a, const T &b)
    {return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}
    
  double pythag(const double a, const double b)
  {
    double absa,absb;

    absa=fabs(a);
    absb=fabs(b);
    if (absa > absb) 
      return absa*sqrt(1.0+sqrt(absb/absa));
    else 
      return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+sqrt(absa/absb)));
  }
  /*
      ** The function
      **                 tqli()
      ** determine eigenvalues and eigenvectors of a real symmetric
      ** tri-diagonal matrix, or a real, symmetric matrix previously
      ** reduced by function tred2[] to tri-diagonal form. On input,
      ** d[] contains the diagonal element and e[] the sub-diagonal
      ** of the tri-diagonal matrix. On output d[] contains the
      ** eigenvalues and  e[] is destroyed. If eigenvectors are
      ** desired z[][] on input contains the identity matrix. If
      ** eigenvectors of a matrix reduced by tred2() are required,
      ** then z[][] on input is the matrix output from tred2().
      ** On output, the k'th column returns the normalized eigenvector
      ** corresponding to d[k]. 
      ** The function is modified from the version in Numerical recipe.
      */
  void tqli(boost::numeric::ublas::vector<double> &d, boost::numeric::ublas::vector<double> &e, boost::numeric::ublas::matrix<double> &z, const int n)
  {
    int m,l,iter,i,k;
    double s,r,p,g,f,dd,c,b;
    
    for (i=1;i<n;i++) e[i-1]=e[i];
    e[n-1]=0.0;
    for (l=0;l<n;l++) {
      iter=0;
      do {
        for (m=l;m<n-1;m++) {
          dd=fabs(d[m])+fabs(d[m+1]);
          if (fabs(e[m])+dd == dd) break;
        }
        if (m != l) {
          iter++;
          try
          {
            iter == 30;
          }
          catch( int e)
          {
            std::cout<<"Too many iterations in tqli"<<std::endl;
          }
          //if (iter++ == 30) nrerror("Too many iterations in tqli");
          g=(d[l+1]-d[l])/(2.0*e[l]);
          r=pythag(g,1.0);
          g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
          s=c=1.0;
          p=0.0;
          for (i=m-1;i>=l;i--) {
            f=s*e[i];
            b=c*e[i];
            e[i+1]=(r=pythag(f,g));
            if (r == 0.0) {
              d[i+1] -= p;
              e[m]=0.0;
              break;
            }
            s=f/r;
            c=g/r;
            g=d[i+1]-p;
            r=(d[i]-g)*s+2.0*c*b;
            d[i+1]=g+(p=s*r);
            g=c*r-b;
            // Next loop can be omitted if eigenvectors not wanted
            for (k=0;k<n;k++) {
              f=z(k,i+1);
              z(k,i+1)=s*z(k,i)+c*f;
              z(k,i)=c*z(k,i)-s*f;
            }
          }
          if (r == 0.0 && i >= l) continue;
          d[l] -= p;
          e[l]=g;
          e[m]=0.0;
        }
      } while (m != l);
    }
  }
}//end NR

using namespace std;
void ReadCommandLine(int argc, char** argv, int &L, float &S)
{
	for (int i=1; i<argc; i++)
	{
		if (strstr(argv[i], "-L" )) {L = atoi(argv[i]+2);};
		if (strstr(argv[i], "-S" )) {S = atof(argv[i]+2);};
	}
}

int main(int argc, char** argv)
{
	int L=12;
	float S=0.5;
  double E_0=1000;
  double E_1=1000;
	double E_ground=2000;
  double E_exc=2000;
  double eps=1e-15;
  int Q=60;//max. dimension of subspace in Q
	int M_S=int(2*S+1);
  const int D=int(pow(M_S,L));//Q Hilbertraum
	ReadCommandLine(argc,argv,L,S);
  
  //basic declarations

  //ublas vector should be known from exercise 1
  boost::numeric::ublas::matrix<double> Vektor(2,D);
  boost::numeric::ublas::vector<double> r(D), u(D);
  
  // We needed those variable, maybe your needs are different?
  int i;
	int index_0;
  double J,JP;
	double Faktor;
	int W=0;
	int X=0;
	int N1,N2,S1,S2,S1P,S2P;
  std::vector<int> P(L);
	//------------------------//
	
  boost::numeric::ublas::vector<double> a(Q),b(Q+1),d(Q),e(Q+1);
	boost::numeric::ublas::matrix<double> z = boost::numeric::ublas::identity_matrix<double>(Q);
  
  // *** Output ****
  std::ofstream outfile;
	std::string fn = str(boost::format("Eigenvalue_spinchain_L=%i_S=%i.txt") % L % S);
	outfile.open(fn);
  
  /* follow this:
   * 1) initialize vectors
   * 2) create random vector r
   * 3) normalize r
   * 4) q_1=r/b_0
   * 5) multiply: H*q_i
   * 6) check boundaries
   * 7) normalize r:
   * 8) solving eigenvalues equation for the tridiagonal system
   *    use: NR::tqli
   * 9) write to file
   *    outfile << W 
              << "\t" 
              << E_0 
              << "\t" 
              << E_1 
              << "\t" 
              << E_1-E_0 
              << endl;
   */
  // initialize vectors

	outfile.close();
}
