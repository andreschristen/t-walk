
/* Program to calculate Integrated Autocorrelation Times: Elisa Reyes-Loza */

#include <math.h>
#include <stdlib.h>		// rand, srand 
#include <stdio.h>	


#ifndef IAT_H
#define IAT_H

#define TINY 1.0e-20


void denom( double x[], double *ax, double *sxx)
{
  unsigned long int j;
  double xt;

  *sxx = 0.0;
  *ax = 0.0;

  for ( j = 0; j < N; j++ )
    *ax += x[j];

  *ax /= (double)N;

  for ( j = 0; j < N; j++ ) {
    xt = x[j] - *ax;
    *sxx += xt * xt;
  }

} /* denom */


double autocorr( double x[], unsigned long int k, double ax, double sxx )
{
  unsigned long int j;
  double sxy = 0.0, r;

  for ( j = 0; j < (N-k); j++ )
    sxy += (x[j] - ax) * (x[j+k] - ax);

  r = sxy / (sxx + TINY);

  return (r);
} /* autocorr */


int IAT( int N, double x[])
{
  unsigned long int i, m = 0;
  double tau = 0, mean, square;
  double sum = 0, Gamma = 0, gamma = 0, gamma_shift = 0;

  denom(x, &mean, &square);

  do { 
    sum += Gamma;
    gamma = autocorr(x, 2*m, mean, square);
    gamma_shift = autocorr(x, 2*m+1, mean, square);
    Gamma = gamma + gamma_shift;
    m ++;
  } while(Gamma > 0);


   return -1 + 2*sum; 
} /* IAT */

#undef TINY

#endif

