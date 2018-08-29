# git
At the introductory to this course you learned how to use git and github.
Here are the basic command again:
```
git clone https://github.com/thierry3000/computer_physics.git
```
# compiling and linking with external libraries
This example is taken from 
[gnu documentation](https://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generator-Examples.html#Random-Number-Generator-Examples)


```
#include <stdio.h>
#include <gsl/gsl_rng.h>

int
main (void)
{
  const gsl_rng_type * T;
  gsl_rng * r;

  int i, n = 10;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for (i = 0; i < n; i++) 
    {
      double u = gsl_rng_uniform (r);
      printf ("%.5f\n", u);
    }

  gsl_rng_free (r);

  return 0;
}
```
We save the this code as `calc_random.cpp`. Trying to compile results in:

>[thierry@dema59 docs]$ gcc -Wall calc_random.cpp -o calc_random
>/tmp/ccd38Bqt.o: In function `main':
>calc_random.cpp:(.text+0x10): undefined reference to `gsl_rng_env_setup'
>calc_random.cpp:(.text+0x17): undefined reference to `gsl_rng_default'
>calc_random.cpp:(.text+0x27): undefined reference to `gsl_rng_alloc'
>calc_random.cpp:(.text+0x46): undefined reference to `gsl_rng_uniform'
>calc_random.cpp:(.text+0x77): undefined reference to `gsl_rng_free'
>collect2: error: ld returned 1 exit status

# python

# exercise 1 - ODE

# exercise 2 - PDE

