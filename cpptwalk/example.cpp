/*
*********

Andres Christen, Ago 2009.

Simple example to use the twalk within a C++ class structure.

This is only useful IF YOU ALREADY KNOW C++.

The twalk info is found at: http://www.cimat.mx/~jac/twalk/

******************** See example.h first please *************************

*
 */


#include <stdio.h>
#include <math.h>

#include "example.h"


int main() {

	// Provide all the parameters:
	
	int dim = 5;
	double  m[5] = { 2.0, 3.9, 5.2, 1.0, 7.8}; //means
	double sd[5] = { 1.2, 2.1, 3.3, 0.1, 2.3}; //std. dev's
	double  a[5] = { 0.0, 0.0, 2.0, 0.0, 3.0}; //Lower boundaries
	
	
	//This sets the seed for for random number generation in the interphase for gsl (in ranfun.h)
	//It most be set before calling TruncIndNor since it uses Unab uniform random number generator.
	Seed((long unsigned int) 0); //to use the time() to get a seed.
	
	//Open a TruncIndNor object
	
	TruncIndNor obj( dim, m, sd, a);
	
	
	//We now open a twalk object:
	
   // Objective function,       initial points,  the obj. fun. dimension, and see below 
	twalk obj_twalk( obj, obj.Getx0(), obj.Getxp0(), obj.get_dim());
	


	//The twalk will use the initial values AND the space provided by x0 and xp0
	//Therefore at the end of the run, we have access to such values.
	
	//And now we run the twalk
	
	obj_twalk.simulation( 100000, (const char *) "twalk.out", "w+",  (int) 1);
	// 100000 iterations, the file to save the output to and if this file should be created etc.
	// (we may use "a+" for example)

	// The last par, save_every, is the iterations saving scheme:
	// if >0, then save iterations that %% save_every == 0 (1 = all iterations)
	// if <0, only *accepted* iterations that %% abs(save_every) == 0 (ie. -1 = all accepeted iter.) 
	// if =0, save all iterations and some info. on each kernel accepatance rates, saved in the file recacc.dat



	
	//The twalk will have dim+1 columns (6 in this case), each for each parameter, and the last one is
	//the corresponding value of eval (-log of the objective).  Very useful for convergence checking.
	//(only the x and not the x' are saved)
	
	//After compiling run example, the output can be seen with R for example, or with CODA etc.
	
	//********************* NOW see the makefile
	
	printf(FAREWELL);
}

