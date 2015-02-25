#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include "auxiliary.h"

#define PI 3.141592653589793

int main (int argc, char* argv[])
{
  int L = pow2(16);
  double dt = 1e-3;
  double t = 1;
  double complex* psi = malloc (L*sizeof(double complex));
  double* V = malloc (L*sizeof(double));
  double* dispersion = malloc (L*sizeof(double));

  fftw_plan psi_forward = fftw_plan_dft_1d (L, psi, 
					    psi, FFTW_FORWARD, FFTW_MEASURE);
  fftw_plan psi_backward = fftw_plan_dft_1d (L, psi, 
					     psi, FFTW_BACKWARD, FFTW_MEASURE);
  for (int i = 0; i < L; ++i)
    {
      psi[i] = 0;
      V[i] = drand ();
      if (i <= L/2)
	dispersion[i] = -2*t*cos(2*PI*i/L);
      else 
	dispersion[i] = -2*t*cos(2*PI*(i-L)/L);
    }
  psi[L/2] = 1;

  for (int i = 0; i < 4000000; ++i)
    {
        fftw_execute (psi_forward);
	for (int j = 0; j < L; ++j)
	  psi[j] *= cexp (-I*dispersion[j]*dt)/L;
	fftw_execute (psi_backward);
	for (int j = 0; j < L; ++j)
	  psi[j] *= cexp (-I*V[j]*dt);
	
	double r2 = 0;
	for (int j = 0; j < L; ++j)
	  r2 += cabs2(psi[j])*square (j-L/2);
	printf ("%g\t%g\n", i*dt, r2);
    }
}

  
  
  
