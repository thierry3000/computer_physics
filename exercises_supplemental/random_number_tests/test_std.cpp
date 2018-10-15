//http://www.cplusplus.com/reference/random/

#include <random>
#include <functional>
#include <iostream>


int main (int argc, char *argv[])
{
  // Seed with a real random value, if available
  std::random_device r;
 
  // use random device as seed
  std::default_random_engine rd(r());
  std::uniform_real_distribution<double> distribution(-1.0, 1.0);
  
  int i, n = 10;
  double value;

  for (i = 0; i < n; i++) 
  {
    std::cout<< distribution(rd) << std::endl;
  }

  //successful return
  return 0;
}
