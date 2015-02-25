#ifndef AUXILIARY_H
#define AUXILIARY_H 1

#include <stdlib.h>
#include <complex.h>
#include <math.h>


/* function square return the square of a real number */
double square (double x);
/* function cabs2 return the square of module of complex number */
double cabs2 (double complex x);
/* function drand generate a uniform random number in [0, 1] */
double drand ();
/* function urand generate a uniform random number in [-1, 1] */
double urand ();
/* function genvec generate a random vector of length S */



/***************************************************/

double square (double x)
{
  return x*x;
}

double cabs2 (double complex x)
{
  return creal (conj(x) *x);
}

double drand ()
{
  return ((double)rand()) / RAND_MAX;
}

double urand ()
{
  return 2 * (drand() - 0.5);
}

int pow2 (int n)
{
  int res = 1;
  while (n-- > 0)
      res <<= 1;
  return res;
}
	
#endif /* AUXILIARY_H */
