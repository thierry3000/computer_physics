/** cip pool compile instructions
 * g++ bosehubbard.cpp -std=c++0x -o run -lboost_program_options -lboost_filesystem -lboost_system
 */

/** some literatur on the topic
 * 
 * @article{BATROUNI199663,
    title = "World line simulations of the bosonic Hubbard model in the ground state",
    journal = "Computer Physics Communications",
    volume = "97",
    number = "1",
    pages = "63 - 81",
    year = "1996",
    note = "High-Performance Computing in Science",
    issn = "0010-4655",
    doi = "http://dx.doi.org/10.1016/0010-4655(96)00022-7",
    url = "http://www.sciencedirect.com/science/article/pii/0010465596000227",
    author = "G.G. Batrouni and R.T. Scalettar",
}
 * @article{PhysRevB.46.9051,
    title = {World-line quantum Monte Carlo algorithm for a one-dimensional Bose model},
    author = {Batrouni, Ghassan George and Scalettar, Richard T.},
    journal = {Phys. Rev. B},
    volume = {46},
    issue = {14},
    pages = {9051--9062},
    numpages = {0},
    year = {1992},
    month = {Oct},
    publisher = {American Physical Society},
    doi = {10.1103/PhysRevB.46.9051},
    url = {https://link.aps.org/doi/10.1103/PhysRevB.46.9051}
}
 *@article{PhysRevLett.65.1765,
    title = {Quantum critical phenomena in one-dimensional Bose systems},
    author = {Batrouni, Ghassan George and Scalettar, Richard T. and Zimanyi, Gergely T.},
    journal = {Phys. Rev. Lett.},
    volume = {65},
    issue = {14},
    pages = {1765--1768},
    numpages = {0},
    year = {1990},
    month = {Oct},
    publisher = {American Physical Society},
    doi = {10.1103/PhysRevLett.65.1765},
    url = {https://link.aps.org/doi/10.1103/PhysRevLett.65.1765}
}
 */


#include "bosehubbard.h"
#define QUICK


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


matrixelements::matrixelements(int N_matrix_elements):N(N_matrix_elements)
{	//constructor
	matrix = new double***[N];
	for(int i = 0; i < N; i++)
  {
  matrix[i] = new double**[N];
	for(int j = 0; j < N; j++)
    {
		matrix[i][j] = new double*[N];
		for(int k = 0; k < N; k++) matrix[i][j][k] = new double[N];
		}
  }
}

matrixelements::~matrixelements()
{			// destruction of the class
	for(int i = 0; i < N; i++)
  {
		for(int j = 0; j < N; j++)
    {
			for(int k = 0; k < N; k++)
      {
				delete [] matrix[i][j][k];
      }
			delete [] matrix[i][j];
    }
		delete [] matrix[i];
  }
	delete [] matrix;
}
/** @brief
 * setting the matrixelements
 */
void matrixelements::setMatrix_weights_of_plaquettes(double hoppingParameter)
{
	for(int i = 0; i < N; i++)  // for matrix elements
		for(int j = 0; j < N; j++)
			for(int k = 0; k < N; k++)
				for(int l = 0; l < N; l++)
        {
					if((i == k) && (j == l)) 
            matrix[i][j][k][l] = exp(- 0.25 * P::deltatau * P::V_0* ( i * (i-1) + j * (j-1) + k * (k-1) + l * (l-1) ));
					else if ((i == k+1) && (l == j+1)) 
            matrix[i][j][k][l] = P::deltatau * hoppingParameter * sqrt(l * (k+1)) * exp(- 0.25 * P::deltatau * P::V_0 * ( i * (i-1) + j * (j-1) + k * (k-1) + l * (l-1) ));
					else if ((i+1 == k) && (l+1 == j)) 
            matrix[i][j][k][l] = P::deltatau * hoppingParameter * sqrt((l+1) * k) * exp(- 0.25 * P::deltatau * P::V_0 * ( i * (i-1) + j * (j-1) + k * (k-1) + l * (l-1) ));
					else matrix[i][j][k][l] = 0.0;
        }
}

void matrixelements::setMatrix_kinetic_energy(double hoppingParameter)
{
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			for(int k = 0; k < N; k++)
				for(int l = 0; l < N; l++)
        {
					if((i == k) && (j == l)) 
            matrix[i][j][k][l] = -1.0 * P::deltatau * hoppingParameter * hoppingParameter * ( k + l + 2.0 * k * l );
					else if ((i == k+1) && (l == j+1)) 
            matrix[i][j][k][l] = - 1.0 / P::deltatau;
					else if ((i+1 == k) && (l+1 == j)) 
            matrix[i][j][k][l] = - 1.0 / P::deltatau;
/*					else if ((i == k+1) && (l == j+1)) matrix[i][j][k][l] = - 1.0 * t * sqrt( (1.0 + k) * l );
					else if ((i+1 == k) && (l+1 == j)) matrix[i][j][k][l] = - 1.0 * t * sqrt( (1.0 + l) * k );
					else if ((i == k+2) && (l == j+2)) matrix[i][j][k][l] = - 1.0 * deltatau * t * t * sqrt( (1.0 + k) * (2.0 + k) * l * (l - 1.0) );
					else if ((i + 2 == k) && (l + 2 == j)) matrix[i][j][k][l] = - 1.0 * deltatau * t * t * sqrt( (1.0 + l) * (2.0 + l) * k * (k - 1.0) );*/
					else 
            matrix[i][j][k][l] = 0.0;
        }
}


/** @brief constructor
 * initialization of the checkerboard and distribution of the bosons on the checkerboard
 */
checkerboard::checkerboard(int Nsite2, int Ntime2, int Nboson2)
{
	Nsite = Nsite2;
	Ntime = Ntime2;
	Nboson = Nboson2;
	int occupancy = Nboson / Nsite;
	int rest = Nboson % Nsite;
	n = new int*[Nsite];
	for(int i = 0; i < Nsite; i++) 
    n[i] = new int[Ntime];
	for(int i = 0; i < rest; i++)
		for(int j = 0; j < Ntime; j++) 
      n[i][j] = occupancy + 1;
	for(int i = rest; i < Nsite; i++)
		for(int j = 0; j < Ntime; j++) 
      n[i][j] = occupancy;
}

checkerboard::~checkerboard(){
	for (int i = 0; i < Nsite ; i++) 
    delete [] n[i];
	delete [] n;
}
/**
 * ####   s    #####
 * we check if a worldline is shifted across is done shaded plaquette, and therefor allowde
 * let (i,j) be the lower left corner of the considered plaquette
 * s = n(i,j) + n(i,j+1) - n(i+1,j) - n(i+1,j+1)
 * and see that
 * s = +2 move from left to right allowed
 * s = -2 move from right to left allowed
 * s not +2 or -2, no move allowed
 * 
 * #### allzero ####
 * check if a worldline hit the considered plaquette at all,
 * if no worldline hits this plaquette there is nothing to do with
 */
void checkerboard::update(double ****matrix){			// local MonteCarlo move, moving a worldline accross an unshaded plaquette
  int i,j, s, move;
	bool sumij, allzero, NotModulo2, notallequal;
  std::array<double,2> final_weights{{0.0,0.0}};
	double p_initial, p_final, zufall = 0.0;
  //this loop finds a valid candidate!
	do 
  {
    double r1 = create_random_number();
    double r2 = create_random_number();
		i = (int) (r1 * Nsite);
		j = (int) (r2 * Ntime);
		sumij = (((i + j) % 2) == 0);
		s = n[(Nsite + i) % Nsite][(Ntime + j) % Ntime] + \
        n[(Nsite + i) % Nsite][(Ntime + j + 1) % Ntime] - \
        n[(Nsite + i + 1) % Nsite][(Ntime + j) % Ntime] - \
        n[(Nsite + i + 1) % Nsite][(Ntime + j + 1) % Ntime];
    
		NotModulo2 = ((abs(s) % 2) == 1);
		allzero = ( (n[(Nsite + i) % Nsite][(Ntime + j) % Ntime] == 0) && 
                (n[(Nsite + i) % Nsite][(Ntime + j + 1) % Ntime] == 0) && 
                (n[(Nsite + i + 1) % Nsite][(Ntime + j) % Ntime] ==0) && 
                ( n[(Nsite + i + 1) % Nsite][(Ntime + j + 1) % Ntime] == 0)
              );
    if (s == 0) 
      notallequal = !(  ( n[(Nsite + i) % Nsite][(Ntime + j) % Ntime] == n[(Nsite + i) % Nsite][(Ntime + j + 1) % Ntime]) && 
                        ( n[(Nsite + i) % Nsite][(Ntime + j + 1) % Ntime] == n[(Nsite + i + 1) % Nsite][(Ntime + j) % Ntime] ) && 
                        (n[(Nsite + i + 1) % Nsite][(Ntime + j) % Ntime] == n[(Nsite + i + 1) % Nsite][(Ntime + j + 1) % Ntime])
                     );
		else 
      notallequal = false;
  }
  while(sumij || allzero || NotModulo2 || notallequal);
  //at this point either of the conditions is fulfilled.
  //and we have a valid candiate at (i,j)
	
	//calculation of the weight of the initial configuration
  //product of matrix elements before move
  p_initial = matrix[n[(Nsite + i - 1) % Nsite][(Ntime + j + 1) % Ntime]]
                    [n[(Nsite + i) % Nsite][(Ntime + j + 1) % Ntime]]
                    [n[(Nsite + i - 1) % Nsite][(Ntime + j) % Ntime]] 
                    [n[(Nsite + i) % Nsite][(Ntime + j) % Ntime]] 
                    * 
              matrix[n[(Nsite + i) % Nsite][(Ntime + j) % Ntime]]
                    [n[(Nsite + i + 1) % Nsite][(Ntime + j) % Ntime]]
                    [n[(Nsite + i) % Nsite][(Ntime + j - 1) % Ntime]]
                    [n[(Nsite + i + 1) % Nsite][(Ntime + j - 1) % Ntime]]
                    *
              matrix[n[(Nsite + i + 1) % Nsite][(Ntime + j + 1) % Ntime]]
                    [n[(Nsite + i + 2) % Nsite][(Ntime + j + 1) % Ntime]] 
                    [n[(Nsite + i + 1) % Nsite][(Ntime + j) % Ntime]] 
                    [n[(Nsite + i + 2) % Nsite][(Ntime + j) % Ntime]] 
                    * 
              matrix[n[(Nsite + i) % Nsite][(Ntime + j + 2) % Ntime]] 
                    [n[(Nsite + i + 1) % Nsite][(Ntime + j + 2) % Ntime]] 
                    [n[(Nsite + i) % Nsite][(Ntime + j + 1) % Ntime]] 
                    [n[(Nsite + i + 1) % Nsite][(Ntime + j + 1) % Ntime]];


	// choosing a move
                    //for the if, only a single move is possible
	if(n[(Nsite + i) % Nsite][(Ntime + j) % Ntime] == 0) 
    move = -1;
	else if (n[(Nsite + i + 1) % Nsite][(Ntime + j) % Ntime] == 0) 
    move = 1;
	else if (n[(Nsite + i) % Nsite][(Ntime + j + 1) % Ntime] == 0) 
    move = -1;
	else if (n[(Nsite + i + 1) % Nsite][(Ntime + j + 1) % Ntime] == 0) 
    move = 1;
	else 
  { //here, two moves are possible, left and right --> check which one is favorable
    int index=0;//required for loop only
    for( move=-1;move<=1;move=move+2)
    {
		final_weights[index] = matrix [n[(Nsite + i - 1) % Nsite][(Ntime + j + 1) % Ntime]] 
                                  [n[(Nsite + i) % Nsite][(Ntime + j + 1) % Ntime] - move] 
                                  [n[(Nsite + i - 1) % Nsite][(Ntime + j) % Ntime]] 
                                  [n[(Nsite + i) % Nsite][(Ntime + j) % Ntime] - move] 
                                  * 
                          matrix [n[(Nsite + i) % Nsite][(Ntime + j) % Ntime] - move] 
                                  [n[(Nsite + i + 1) % Nsite][(Ntime + j) % Ntime] + move] 
                                  [n[(Nsite + i) % Nsite][(Ntime + j - 1) % Ntime]] 
                                  [n[(Nsite + i + 1) % Nsite][(Ntime + j - 1) % Ntime]] 
                                  *
                          matrix [n[(Nsite + i + 1) % Nsite][(Ntime + j + 1) % Ntime] + move] 
                                  [n[(Nsite + i + 2) % Nsite][(Ntime + j + 1) % Ntime]] 
                                  [n[(Nsite + i + 1) % Nsite][(Ntime + j) % Ntime] + move] 
                                  [n[(Nsite + i + 2) % Nsite][(Ntime + j) % Ntime]] 
                                  * 
                          matrix [n[(Nsite + i) % Nsite][(Ntime + j + 2) % Ntime]] 
                                  [n[(Nsite + i + 1) % Nsite][(Ntime + j + 2) % Ntime]] 
                                  [n[(Nsite + i) % Nsite][(Ntime + j + 1) % Ntime] - move] 
                                  [n[(Nsite + i + 1) % Nsite][(Ntime + j + 1) % Ntime] + move];
      index=index+1;
    }
    if(final_weights[1]>final_weights[0])
      move = 1;
		else 
      move = -1;
  }

	// calculation of the weight of the final configuration
  p_final = matrix  [n[(Nsite + i - 1) % Nsite][(Ntime + j + 1) % Ntime]] 
                    [n[(Nsite + i) % Nsite][(Ntime + j + 1) % Ntime] - move] 
                    [n[(Nsite + i - 1) % Nsite][(Ntime + j) % Ntime]] 
                    [n[(Nsite + i) % Nsite][(Ntime + j) % Ntime] - move] 
                    * 
            matrix  [n[(Nsite + i) % Nsite][(Ntime + j) % Ntime] - move] 
                    [n[(Nsite + i + 1) % Nsite][(Ntime + j) % Ntime] + move] 
                    [n[(Nsite + i) % Nsite][(Ntime + j - 1) % Ntime]] 
                    [n[(Nsite + i + 1) % Nsite][(Ntime + j - 1) % Ntime]] 
                    * 
            matrix  [n[(Nsite + i + 1) % Nsite][(Ntime + j + 1) % Ntime] + move] 
                    [n[(Nsite + i + 2) % Nsite][(Ntime + j + 1) % Ntime]] 
                    [n[(Nsite + i + 1) % Nsite][(Ntime + j) % Ntime] + move] 
                    [n[(Nsite + i + 2) % Nsite][(Ntime + j) % Ntime]] 
                    * 
            matrix  [n[(Nsite + i) % Nsite][(Ntime + j + 2) % Ntime]] 
                    [n[(Nsite + i + 1) % Nsite][(Ntime + j + 2) % Ntime]] 
                    [n[(Nsite + i) % Nsite][(Ntime + j + 1) % Ntime] - move] 
                    [n[(Nsite + i + 1) % Nsite][(Ntime + j + 1) % Ntime] + move];

	// drawing a random number
	zufall = create_random_number();

	// Metropolis step, update of the configuration
	if( zufall < p_final/p_initial)
  {
		n[(Nsite + i) % Nsite][(Ntime + j) % Ntime] -= move;
		n[(Nsite + i) % Nsite][(Ntime + j + 1) % Ntime] -= move;
		n[(Nsite + i + 1) % Nsite][(Ntime + j) % Ntime] += move;
		n[(Nsite + i + 1) % Nsite][(Ntime + j + 1) % Ntime] += move;
  }
}

void checkerboard::printworldlines(std::ofstream &daten)
{  
	for(int j = 0; j < Ntime; j++){
		for(int i = 0; i < Nsite; i++) daten << n[i][j] << "\t";
		daten << std::endl;
		}
	daten << std::endl;
}

/** @brief
 * constructor
 */
correlationtime::correlationtime(checkerboard &worldlines, double ****matrix, int tmax2):tmax(tmax2)
{
  int Nsite = worldlines.getNsite();
  int Ntime = worldlines.getNtime();
  
  ///         INIT
	n = new int**[tmax];
	for(int i = 0; i < tmax; i++)
  {
		n[i] = new int *[Nsite];
		for(int j = 0; j < Nsite; j++)
    {
			n[i][j] = new int[Ntime];
			for(int k = 0; k < Ntime; k++) 
        n[i][j][k] = 0;
    }
  }
	n2 = new double**[tmax/2];
	for(int i = 0; i < tmax/2; i++)
  {
		n2[i] = new double *[Nsite];  
		for(int j = 0; j < Nsite; j++)
    {
			n2[i][j] = new double[Ntime];
			for(int k = 0; k < Ntime; k++) n2[i][j][k] = 0.0;
		}
	}
	for(int j = 0; j < Nsite; j++)
		for(int k = 0; k < Ntime; k++) n2[0][j][k] = 1.0;
	m = new double[tmax/2];
	for(int i = 0; i < tmax/2; i++) m[i] = 0.0;
  
  // initialize the array n
  for(int i = 0; i < tmax; i++)        //for all timepoints we consider
  {
    worldlines.update(matrix);
    for(int j = 0; j < Nsite; j++)     //for every site
      for(int k = 0; k < Ntime; k++)   // for every trotter time point
        n[i][j][k] = worldlines.n[j][k];
  }
}

correlationtime::~correlationtime(){		//destruction of the class
	for(int i = 0; i < tmax; i++){
		for(int j = 0; j < P::Nsites; j++) delete [] n[i][j];
		delete [] n[i];
		}
	delete [] n;

	for(int i = 0; i < tmax/2; i++){
		for(int j = 0; j < P::Nsites; j++) delete [] n2[i][j];
		delete [] n2[i];
		}
	delete [] n2;
	delete [] m;
}

/** @brief
 * prints the correlations as a timedependent function
 */
void correlationtime::print(std::ofstream &correlation)
{
	double m1, m2, m3;
	for(int zeit = 0; zeit < tmax/2; zeit++)  // for all zeit we consider
  {
		if ((zeit % 25) == 0)                   // calculate correlation every 25 steps
    {
		for(int i = 0; i < P::Nsites; i++)          //for every site
			for(int j = 0; j < P::Ntrotter; j++)        //for every trotter time point
      {
				m1 = 0.0;
				m2 = 0.0;
				m3 = 0.0;
				for(int t2 = 0; t2 < tmax - zeit; t2++)  // we only need to sum until we there
        {
					m1 += 1.0 * n[t2][i][j] * n[t2+zeit][i][j];
					m2 += 1.0 * n[t2][i][j];
					m3 += 1.0 * n[t2+zeit][i][j];
				}
				m1 = 1.0 * m1 / (tmax-zeit);
				m2 = 1.0 * m2 / (tmax-zeit);
				m3 = 1.0 * m3 / (tmax-zeit);
				n2[zeit][i][j] = (m1 - m2 * m3); 
      }
		for(int i = 0; i < P::Nsites; i++)
			for(int j = 0; j < P::Ntrotter; j++) 
        m[zeit] += n2[zeit][i][j];
		m[zeit] = m[zeit] / (1.0 * P::Nsites) / (1.0 * P::Ntrotter); // normalize
		correlation << zeit << "\t" << m[zeit] / m[0] << std::endl;
		}
  }
}

/** @brief 
 * calculate kinectic energy
 */
double kinenergy(double ****kinetic1, int **n, int Nsite, int Ntime)
{
	double kineticenergy = 0.0;

	for(int j = 0; j < Ntime; j++)
  {
		for(int i = 0; i < Nsite; i++)
    {
			if( (i%2) == 0 ) 
        kineticenergy += kinetic1 [n[i][(j+1)%Ntime]]
                                  [n[(i+1)%Nsite][(j+1)%Ntime]]
                                  [n[i][j]][n[(i+1)%Nsite][j]] ;
			else 
        kineticenergy += kinetic1 [n[i][j]]
                                  [n[(i+1)%Nsite][j]]
                                  [n[i][(Ntime+j-1)%Ntime]]
                                  [n[(i+1)%Nsite][(Ntime+j-1)%Ntime]];
		}
	}
	kineticenergy = kineticenergy / Ntime;
	return kineticenergy;
}

kinetic::kinetic(int Nsite2, int Ntime2, int Nboson2, int correltime2):Nsite(Nsite2), Ntime(Ntime2), Nboson(Nboson2), correltime(correltime2)
{
}

void kinetic::print(double ****matrix){
  double kin_energy= calculate(matrix);
  std::cout << "bosondensity   " << 1.0 * Nboson / Nsite <<  "    kinetic energy    " << kin_energy << std::endl;
}

double kinetic::calculate(double ****matrix)
{
	matrixelements kinetic1(25);
	kinetic1.setMatrix_kinetic_energy(P::t);
	double kineticenergy;
	checkerboard worldlines(Nsite, Ntime, Nboson);
  int iters =0;
#ifdef QUICK
  iters = 1000;
#else
  iters = 10000;
#endif
	for(int i = 0; i < iters; i++) 
    worldlines.update(matrix);
	kineticenergy = kinenergy(kinetic1.matrix, worldlines.n, Nsite, Ntime);
  return kineticenergy;
}


int main(int argc, char **argv)
{
#ifdef OLD_RANDOM
  std::cout << "OLD_RANDOM defined " << std::endl;
#endif
  handle_program_args(argc, argv);
  std::string outputFolderName = "QMC_output";
  boost::filesystem::path dir(outputFolderName);
  if(!(boost::filesystem::exists(dir)))
  {
    std::cout<<"QMC_output folder does not exist."<<std::endl;
    if (boost::filesystem::create_directory(dir))
      std::cout << "....Successfully Created !" << std::endl;
  }
  std::cout << "creating matrix elements! " << std::endl;
	matrixelements H1(P::NMatrixElements);
  std::cout << " inital set matrix " << std::endl;
	H1.setMatrix_weights_of_plaquettes(P::t);
  std::cout << "matrix set " << std::endl;
  
  checkerboard bosonworldlines(P::Nsites, P::Ntrotter, P::Nboson);
  std::cout << " checkerboard initialized " << std::endl;
  
  std::ofstream worldlines(outputFolderName + "/" + P::fn); // the output will to there
  bosonworldlines.printworldlines(worldlines);
  
  std::cout << " printed initial bosonworldlines " << std::endl;
#ifdef QUICK
  std::cout << "starting quick run " << std::endl;
  for(int i = 0; i < 1000; i++)
    bosonworldlines.update(H1.matrix);
#else
  for(int i = 0; i < 10000; i++)
    bosonworldlines.update(H1.matrix);
#endif
  bosonworldlines.printworldlines(worldlines);
  std::cout << " printed altered bosonworldlines " << std::endl;
  worldlines.close();


	return 0;
}
