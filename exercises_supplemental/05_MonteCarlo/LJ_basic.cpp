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
 * g++ LJ_basic.cpp -std=c++0x -o myBarometricFormula -lboost_program_options
 */

#include<string.h>
#include<math.h>
#include<iostream>
#include<fstream>
#include <vector>
#include<list>

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

using namespace std;

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

struct vektor2d
{
double x;
double y;
};

typedef std::vector<vektor2d> zustand;

//model parameters (default)
namespace P // Parameters
{
  //model
double beta = 10.01,
       m = 0.01,
       eps = 1.0,
       sigma = 0.2,
       a = 10.0,
       b = 10.0;
int    N = 30;
  //simulation
long int steps=pow(10,4);
int    burnIn=1000;
}

//general program options
string output_filename = "Lenard_Jones";
static double dxmax = 0.1*P::a;
static double dymax = 0.1*P::b;
int movecounter=0;

void print_description(std::ostream &os)
{
  os << "Simulation of Lennard- Jones Potential with gravity:" << endl;
  os << "  H = \\sum m*g*x_i + " << endl;
  os << "  \\sum_i \\sum_{j=i+1} 4\\epsion((\\sigma/(x_i-x_j))^12-(\\sigma/(x_i-x_j))^6)" << endl;
  os << "  Computerphysics lecture 2017, AG Rieger" << endl;
}


/* why we do this?????
example:
compile with:
g++ LJ_basic.cpp -std=c++0x -O3 -o run -lboost_program_options

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
    ADD_PARAM("beta", P::beta, "Inverse temperature")
    ADD_PARAM("m", P::m, "mass of particles")
    ADD_PARAM("eps", P::eps, "Lennard- Jones parameter")
    ADD_PARAM("sigma", P::sigma, "Lennard- Jones parameter")
    ADD_PARAM("a", P::a, "extension first dimension")
    ADD_PARAM("b", P::b, "extension second dimension")
    ADD_PARAM("N", P::N, "number of particles")
    ADD_PARAM("steps", P::steps, "number of Monte Carlo steps")
    ADD_PARAM("burnIn", P::burnIn, "neglected number of steps")
    ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (vm.count("help"))
  {
    cout << desc << endl;
    cout << "The following is done here:" << endl;
    print_description(cout);
    exit(1);
  }
  cout << "Monte Carlo run "; 
  print_description(cout);
  cout << "Model Parameters:" << endl;
  #define OUT_PARAM1(PARAM) cout << #PARAM" = " << P::PARAM << endl
  OUT_PARAM1(N);
  OUT_PARAM1(beta);
  OUT_PARAM1(steps);
  OUT_PARAM1(sigma);
  OUT_PARAM1(eps);
  OUT_PARAM1(a);
  OUT_PARAM1(b);
  OUT_PARAM1(m);
  OUT_PARAM1(burnIn);
  
  cout << "output-filename: " << output_filename << endl;
}

/*              output to file (c)omma (s)eperated (v)alue file

This is helpfull if you like to analyze your data with a more sophisticated data analyzing tools like
MATLAB, Mathematica or freemat
First column is the x coordinate, second y coordinate.

Of course your are free to use anything you like to create nice histograms or plot ;-)
*/

void writeToCSV(list<zustand> &pzustandsliste)
{
  string fn = (boost::format("%s_N_%i_steps_%i.csv")
                  % output_filename
                  % P::N
                  % P::steps
              ).str();
  std::ofstream f(fn.c_str());
  if (!f.good())
    throw std::runtime_error("unable to open output file: "+fn);

  // see http://www.cplusplus.com/reference/list/list/
  std::list<zustand>::const_iterator iterator_over_zustaende;
  std::vector<vektor2d>::const_iterator iterator_over_positions;

  for(iterator_over_zustaende = pzustandsliste.begin(); iterator_over_zustaende != pzustandsliste.end(); ++iterator_over_zustaende)
  {//loop over states
    for(iterator_over_positions =(*iterator_over_zustaende).begin();
      iterator_over_positions != (*iterator_over_zustaende).end();
      ++iterator_over_positions)//loop over positions
    {
      f << (*iterator_over_positions).x << "," << (*iterator_over_positions).y << std::endl;
    };
  };
};//end writeToCSV


double d2(const vektor2d & v1, const vektor2d & v2) //distance^2
{
	//both left
	vektor2d delta;
	delta.x= v1.x-v2.x;
	delta.y= v1.y-v2.y;
	return sqrt(pow(delta.x,2)+pow(delta.y,2));
};

//Lennard Jones
double V(const vektor2d & v1, const vektor2d & v2) 
{
  return 4*P::eps*( pow(pow(P::sigma,2)/d2(v1,v2),6)-pow(pow(P::sigma,2)/d2(v1,v2),3)); //Minimum bei 2^(1/6)*sigma
};

double V_op (const vektor2d & v1)
{
   return P::m*9.81*v1.y;
}

//total energy
double H(const zustand & z) 
{
  int i,j;
  double sum=0.0;
  for(i=0;i<z.size()-1;++i)
    {
    for(j=i+1;j<z.size();++j)
      {
      sum+=V(z[i],z[j]);
      };
    };
  return sum;
};

//total energy
double Hp_op(const zustand & z) 
{
int i,j;
double sum=0.0;
for(i=0;i<z.size()-1;++i)
        {
         sum+=V_op(z[i]);
	/*
         for(j=i+1;j<z.size();++j)
                {
                 sum+=V(z[i],z[j]);
                };
	*/
        };
return sum;
};



zustand Zustandwuerfeln(zustand z,int & nr)
{
  double dx;
  double dy;
  int move_nr=floor((create_random_number())*P::N);
  if (move_nr==P::N) --move_nr;
  nr=move_nr;

  dx=(2.0*(create_random_number())-1.0)*dxmax;
  dy=(2.0*(create_random_number())-1.0)*dymax;

  vektor2d buffer;
  buffer.x=z[nr].x+dx;
  buffer.y=z[nr].y+dy;

  //non periodic
  if((buffer.x<0||buffer.x>P::a)||(buffer.y<0||buffer.y>P::b))
    {
    buffer.x=z[nr].x;
    buffer.y=z[nr].y;
    z[nr]=buffer;
    }
    
  z[nr]=buffer;

return z;


//periodic
/*should be implemented in exercise 1b) */

}

/* exercise 1c) */
void density()
{
  
}

zustand Initialisierung()
{
zustand z(P::N);
for(int i=0;i<P::N;++i)
{
	z[i].x=(create_random_number())*P::a;
	z[i].y=(create_random_number())*P::b;	
};
return z;
}

zustand Metropolis_step(zustand z)
{
int nr;
zustand z_neu=Zustandwuerfeln(z, nr);
if (exp(-P::beta*(Hp_op(z_neu)-Hp_op(z)))>=(create_random_number())) {++movecounter; return z_neu; }
//exercise 1 a)
//if (exp(-beta*dE(z_neu,z,nr))>=(create_random_number())) {++movecounter; return z_neu; }

else return z;
}




int main(int argc, char** argv) 
{
  handle_program_args(argc, argv);
  //Initialisierung
  list<zustand> zustandsliste;
  zustand pos=Initialisierung();

  //Monte-Carlo
  for(unsigned long i=1;i<=P::steps;++i)
  {
    if(i>P::burnIn)
    {
      zustandsliste.push_back(pos);
    }
    pos=Metropolis_step(pos);	
  };
  writeToCSV(zustandsliste);
}  // end main 
