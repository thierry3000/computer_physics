# linux OS
## basic commands
We assume that you are familiar with basic linux command such as
```
ls, cd, rm, .., echo, chmod, mkdir, find, etc
```
Never the less, let us

- [x] create a folder name "computer_physics_2018"
- [x] create a text file with your name in that folder e.g. "Thierry_Fredrich.txt" in that folder
- [x] write your student id (Matrikelnummer), eMail address and type of study into the file
- [x] (secure) copy this file to the main cip server (gordon) to the folder /home/comphys/students_2018
- [x] delete the folder and the text file

To obtain the latest source, you can either visit the [github page](https://github.com/thierry3000/computer_physics) and hit the download button or you use a git tool to get your own working copy of the repository. We do not have the time to explain all of it, but the following line will copy the code.

```
git clone https://github.com/thierry3000/computer_physics.git
```
## Bash scripting (skip until spare time)
Bash is one of the languages that allows direct communication with the linux kernel that is the heart of a linux based OS. Note that there are other languages as well. All commands typed into the terminal could be saved in a text file which is than called bash script. We do a simple example. Run the following line in your console

```
for i in {0..10}; do echo ${i}; done
```
than save the line in a text file called "loop.sh". By convention we have to add and additional line in the beginning:
```
#!/bin/bash
```
the shebang tells the machine the language of your commands. In this case we speak "bash".

- [x] Allow execution of your text file.
- [x] Execute your script.

# compiling and linking with external libraries

## basic example
This example is taken from [gnu documentation](https://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generator-Examples.html#Random-Number-Generator-Examples)


```
#include <stdio.h>        //printf
#include <gsl/gsl_rng.h>  //gsl random

/**
 * main is the interface from your code to 
 * the operating system!!!
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
```
We save the this code as `calc_random.cpp`. Trying to compile with the command 

>[thierry@dema59 docs]$ g++ calc_random.cpp -Wall -o calc_random

results in:

>/tmp/ccd38Bqt.o: In function main:
>calc_random.cpp:(.text+0x10): undefined reference to `gsl_rng_env_setup'
>calc_random.cpp:(.text+0x17): undefined reference to `gsl_rng_default'
>calc_random.cpp:(.text+0x27): undefined reference to `gsl_rng_alloc'
>calc_random.cpp:(.text+0x46): undefined reference to `gsl_rng_uniform'
>calc_random.cpp:(.text+0x77): undefined reference to `gsl_rng_free'
>collect2: error: ld returned 1 exit status

which means that we forgot to provide the implementation for gsl. We 
correct to:

>[thierry@dema59 docs]$ g++ calc_random.cpp -Wall -lgsl -o calc_random

which gives an even longer list of error. This time the errors look like

> /usr/lib/gcc/x86_64-pc-linux-gnu/8.1.0/../../../../lib/libgsl.so: undefined reference to `cblas_ztrsv'

letting you know that there is an other undefined reference to something with cblas. 
We add the link to cblas

>[thierry@dema59 docs]$ g++ calc_random.cpp -Wall -lgsl -lgslcblas -o calc_random

and the compilation proceeds. Now we can execute our binary by

> ./calc_random

## include
2 lines of include are stated in the example above.

```
#include <stdio.h>
#include <gsl/gsl_rng.h>
```
This tells the compiler to copy the content of those to files into our example. But how does the compiler know where to look for those files?
Dependent on you operating system, some default locations are scanned. For most linux distributions those locations are:

- /include
- /usr/include
- /usr/local/include

The rectangular brackets (<>) are reminiscent of using the default locations. 
Apart from the default locations you can always add non-default locations 
with the option "-I" to the compiler or provide a relative or absolute path on your file system. 

- [x] Try to find the files "stdio.h" and "gsl_rng.h" on your system!

## linking
Linking is different to include and therefore done at a different step --- the linking procedure.
In contrast to "include" which just copies the content of a file,
the linking stage guides your program to locations where other 
machine usable code could be found and directly used.

In our compile line we added "-lgsl" and "-lgslcblas".
Akin to include, the operating system provides default locations for 
libraries which is how we call the directly usable machine code.

- /lib
- /usr/lib 
- /usr/lib64
- /usr/lib/x86_64-linux-gnu (different architectures might we present on the same file system) 
- /usr/local/lib 
- /usr/local/lib64

On unix based systems, shared libraries follow the naming convention that 
they start with "lib" and end with ".so" (we will not consider other types here). On an other popular 
operating system they are call ".dll" files.

- [x] Try to find the files "libgsl.so" and "libgslcblas.so" on your system!
- [x] Try to compile and run the main function within the "hello_world.cpp"

# optimization and debugging
Yet an other option of the compiler is related to the level of optimization. Without any additional flags the optimization level is "-O0". To speed up your code you can increase this level up to "-O3". Applied to our recent example this would result in the command

> g++ calc_random.cpp -Wall -lgsl -lgslcblas -o calc_random -O3

- [x] If you have other code at hand, compare the runtimes with and without the "-O3" flag.

Optimization is one end of the ladder, but there is an other direction, too. Assume you are hunting for a runtime error and do not care so much about optimization, but you are more interested in a correct execution. Then you should tell the compiler to add additional debugging information with "-g" flag.

Assume, for some reason, that we did something wrong in our hello world example (uncomment the line) and compile again with "-g" flag. Run the gnu debugger gdb with the name of your binary as argument:

> gdb your-binary-name

type "run" at hit enter.


See [here](https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html) and [here](https://gcc.gnu.org/onlinedocs/gcc/Debugging-Options.html#Debugging-Options) for all compiler options.

# python

In contrast to compile the sources to specific hardware, there is an other approach called interpreter based programming. This is widely used  in science because it is very practical (Matlab and friends). 

- [x] What are the pros and cons of compiler based vs. interpreter based approaches?

## interpreter

Open a console and type:

> python

This will open an interpreter and you can directly enter your commands.

## scripting

An other handy way of using python is to write scripts and hand them to the interpreter.

> python hello_world.py

- [x] Try to understand the source of hello_world.py
- [x] Plot a sine.

## packages

Please have a look at the following folder:

> /usr/lib/python2.7/dist-packages

- [x] What can you find there?
- [x] What is a __init__.py file?

# IDE (if time left)

# exercise 1 - ODE

# exercise 2 - PDE

