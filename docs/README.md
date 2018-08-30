# git
At the introductory to this course you learned how to use git and github.
Here are the basic command again:
```
git clone https://github.com/thierry3000/computer_physics.git
```
# compiling and linking with external libraries

## basic example
This example is taken from [gnu documentation](https://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generator-Examples.html#Random-Number-Generator-Examples)


```
#include <stdio.h>
#include <gsl/gsl_rng.h>

int main (void)
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
We save the this code as `calc_random.cpp`. Trying to compile with the command 

>[thierry@dema59 docs]$ g++ calc_random.cpp -Wall -o calc_random

results in:

>/tmp/ccd38Bqt.o: In function `` ` ``main':
>calc_random.cpp:(.text+0x10): undefined reference to `` ` ``gsl_rng_env_setup'
>calc_random.cpp:(.text+0x17): undefined reference to `` ` ``gsl_rng_default'
>calc_random.cpp:(.text+0x27): undefined reference to `` ` ``gsl_rng_alloc'
>calc_random.cpp:(.text+0x46): undefined reference to `` ` ``gsl_rng_uniform'
>calc_random.cpp:(.text+0x77): undefined reference to `` ` ``gsl_rng_free'
>collect2: error: ld returned 1 exit status

which means that we forgot to provide the implementation for gsl. We 
correct to:

>[thierry@dema59 docs]$ g++ calc_random.cpp -Wall -lgsl -o calc_random

which gives an even longer list of error. This time they look like

> /usr/lib/gcc/x86_64-pc-linux-gnu/8.1.0/../../../../lib/libgsl.so: undefined reference to `cblas_ztrsv'

letting you know that there is an other undefinde reference to something with cblas. 
We add the link to cblas

>[thierry@dema59 docs]$ g++ calc_random.cpp -Wall -lgsl -lcblas -o calc_random

and the compilation proceeds.

## include
2 lines of include are stated in the code above.

```
#include <stdio.h>
#include <gsl/gsl_rng.h>
```
This tells the compiler to copy the content of those to files into our example. But how does the compiler know where to look for thos files?

## linking


# python

# exercise 1 - ODE

# exercise 2 - PDE

