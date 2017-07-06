#ifndef BOSEHUBBARD_H
#define BOSEHUBBARD_H  // header guard

#include<iostream>   // zum Einbinden der Standard Ein und Ausgabe
#include<math.h>     // beiinhaltet die Standard - Mathe Funktionen ( Wurzel, sin, cos ... )
#include<fstream>    // zum Schreiben / Lesen aus in Dateien
#include<iomanip>    // Zum Arbeiten mit Strings
#include<cstdlib>    // beiinhaltet Zufallsgenerator
#include<time.h>
#include<cmath>

// c++ version von printf.
#include <boost/format.hpp>
// Boost program_options support command line arguments.
#include <boost/program_options.hpp>
// allows us to create folders
#include <boost/filesystem.hpp>

#include<gsl/gsl_errno.h> // error handling of gsl --> maybe you get some information for free
#include<gsl/gsl_fft_real.h> // computes fourier transforms
// handels smart pointers like std::shared_ptr
#include <memory>

//model parameters (default)
namespace P // Parameters
{
  //model
double deltatau = 0.0625,
       t = 1.0,
       V_0 = 20.0;
int    Nsites = 16,
       Ntrotter = 64,
       Nboson = 16,
       NMatrixElements = 25;
std::string fn = "Worldlines.dat";
}

void print_description(std::ostream &os)
{
  os << "Simulation of Bose hubbard model" << std::endl;
  os << "       ground state, periodic boundary conditions" <<std::endl;
  os << "  Computerphysics lecture 2017, AG Rieger" << std::endl;
}

/* why we do this?????
example:
compile with:
g++ LJ_basic.cpp -O3 -o run -lboost_program_options

then you can set the parameter during runtime using, see also exercise 2, PDE

./run --help             shows options

*/
void handle_program_args(int argc, char **argv)
{  
  namespace po = boost::program_options;

  // Makro für Program Optionen in die options_description Liste zu machen.
  #define ADD_PARAM(NAME, P, MSG) (NAME, po::value(&P)->default_value(P), MSG)
  po::options_description desc("Options");
  desc.add_options()
    ("help", "generate help message")
    ADD_PARAM("deltatau", P::deltatau, "imaginary timestep")
    ADD_PARAM("Nsites", P::Nsites, "number of sites")
    ADD_PARAM("Ntrotter", P::Ntrotter, "number of trotterslices")
    ADD_PARAM("Nboson", P::Nboson, "number of bosons")
    ADD_PARAM("t", P::t, "Hopping parameter")
    ADD_PARAM("V0", P::V_0, "Onsite potential")
    ADD_PARAM("fn", P::fn, "output file name")
    ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (vm.count("help"))
  {
    std::cout << desc << std::endl;
    std::cout << "The following is done here:" << std::endl;
    print_description(std::cout);
    exit(1);
  }
  print_description(std::cout);
  std::cout << "Model Parameters:" << std::endl;
  #define OUT_PARAM1(PARAM) std::cout << #PARAM" = " << P::PARAM << std::endl
  OUT_PARAM1(deltatau);
  OUT_PARAM1(Nsites);
  OUT_PARAM1(Ntrotter);
  OUT_PARAM1(Nboson);
  OUT_PARAM1(t);
  OUT_PARAM1(V_0);
  //OUT_PARAM1(fn);
  
  std::cout << "output-filename: " << P::fn << std::endl;
}

// declare stuff needed for bosehubbard
double kinenergy(double ****kinetic1, int **n, int Nsite, int Ntime);

/** @brief
 * class for checkerboard decomposition of the Hamiltonian
 */
class checkerboard
{
  public:
    int **n;
    checkerboard(int Nsite2, int Ntime2, int Nboson);
    ~checkerboard();
    void update(double ****matrix);
    void printworldlines(std::ofstream &); // streaming the worldlines to a file
    int getNsite(){return Nsite;}
    int getNtime(){return Ntime;}
    int getNboson(){return Nboson;}
  private:
    int Nsite;		// number of sites
    int Ntime;		// number of timeslices
    int Nboson;		// number of bosons
};

// /** @brief 
//  * calculation of kinetic energy terms
//  */
// class matrixelements2{
// public:
// 	matrixelements2(int N_total_boson);
// 	~matrixelements2();
// 	void setMatrix();
// 	double ****matrix;
// private:
// 	int N;
// };

/** @brief
 * Class to calculate the weights of the plaquettes
 */
class matrixelements{
public:
	matrixelements(int N_matrix_elements);
	~matrixelements();
	void setMatrix_weights_of_plaquettes(double hoppingParameter);
  void setMatrix_kinetic_energy(double hoppingParameter);
	double ****matrix;
private:
	int N; // number of matrix elements
};

/** @brief 
 * calculation of correlations of MonteCarlo Configurations
 */
class correlationtime
{
private:
  int tmax;
  int ***n;
  double ***n2;
  double *m;

public:
  correlationtime(checkerboard&, double ****matrix, int tmax2);
  ~correlationtime();
  void initialize(checkerboard &, double ****matrix);
  void print(std::ofstream &correlation);
};

class kinetic{
private:
  int Nsite;
  int Ntime;
  int Nboson;
  int correltime;

public:
	kinetic(int Nsite2, int Ntime2, int Nboson2, int correltime2);
  double calculate(double ****matrix);
	void print(double ****matrix);
  int getNsite(){return Nsite;}
  int getNtime(){return Ntime;}
  int getNboson(){return Nboson;}
  int getCorreltime(){return correltime;}
};

/// solution
namespace solution
{
  class density{
  private:
    int Nsite;
    int Ntime;
    int Nboson;
    int correltime;

  public:
    density(int correltime2, int Nsite2, int Ntime2);
    void print(double ****matrix, double &t, double &deltatau, std::ofstream &);
  };
  
  class phasediagram
  {
  private:
    int Nsite;
    int Ntime;
    int Nboson;
    int correltime;

  public:
    phasediagram(int correltime2, int Nsite2, int Ntime2);
    void print(double ****matrix, double &t, double &deltatau, std::ofstream &);
  };
  
  /** @brief 
  * class for the calculation of the pseudo current current correlation
  */
  class currentcorrelation
  {
  private:
    const int correltime;
    const int Nsite;
    const int Ntime;
    std::vector<double> *J;
    std::vector<double> *J2;

  public:
    currentcorrelation(const int correltime2, const int Nsite2, const int Ntime2);
    void calculate(checkerboard &, double ****matrix);
    const int getNtime() const {return Ntime;};
  };
}

#endif
