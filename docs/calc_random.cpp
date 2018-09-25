#include <stdio.h>        //printf
#include <gsl/gsl_rng.h>  //gsl random

/**
 * main is the interface from your code to 
 * the operating sytem!!!
 * 
 * by convention, it returns an integer as 
 * return state.
 * successful exit: 0
 */
int main (void)
{
  //reserve space for pointers
  const gsl_rng_type * T;
  gsl_rng * r;

  int i, n = 10;
  
  //standard initialization
  gsl_rng_env_setup();

  //allocate pointers
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for (i = 0; i < n; i++) 
    {
      //returns a random number
      double u = gsl_rng_uniform (r);
      //prints float with 5 digets
      printf ("%.5f\n", u);
    }

  gsl_rng_free (r);
  //successful return
  return 0;
}
