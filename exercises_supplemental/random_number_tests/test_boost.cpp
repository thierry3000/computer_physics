#include <boost/random.hpp>
#include <boost/random/random_device.hpp>

boost::random_device rd;
boost::random::mt19937_64 gen(rd());
//https://www.boost.org/doc/libs/1_66_0/doc/html/boost_random/reference.html#boost_random.reference.distributions
boost::random::normal_distribution<double> dis;
double value = dis(gen);

int main (void)
{
  int i, n = 10;
  double value;

  for (i = 0; i < n; i++) 
    {
      //returns a random number
      double u = dis(gen);
      //prints float with 5 digets
      printf ("%.5f\n", u);
    }

  //successful return
  return 0;
}
