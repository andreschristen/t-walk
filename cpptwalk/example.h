/*
*********

Andres Christen, Oct 2012.

Simple example to use the twalk within a C++ class structure.

This is only useful IF YOU ALREADY KNOW C++.

The twalk info is found at: http://www.cimat.mx/~jac/twalk/

**********/



#ifndef EXAMPLE_H
#define EXAMPLE_H 

//This is my traditional farewell, may be changed to something more "serious" or to ""
#define FAREWELL "Eso es to...eso es to...eso es to...eso es toooodo amigos!\n"

#include <stdio.h>
#include <math.h>

#include "ranfun.h"
#include "twalk.h"              // twalk simulator class

/*****************************
We are going to run the twalk with a simple objective function:
n independent normal distributions, but each truncated from below at a[j]
*/


/*****************************
We create a class that derives from the templete class obj_fcn (objective function) defined in twalk.h
and this is the interphase with the twalk.  Some few functions are needed to be define
to create an instance of obj_fcn .
*/

class TruncIndNor: public obj_fcn {

	private:
		
		double *m; //pointer to the set of mean's for each normal
		double *sd; //pointer to the array of standard dev's for each normal
		double *a; //the lower boundaries for the support of each normal
		
		double *x0, *xp0; //two initial points

	public:
		//Our constructor.  The obj_fcn constructor takes one int
		//as argument, dim, which is the dimension of the objective function 
		//we take that as well an pass it on.
		TruncIndNor( int dim, double *mm, double *ssd, double *aa) : obj_fcn(dim) {
		
			//This is the rest of our constructor
			
			//Create space
			m  = new double[dim];
			sd = new double[dim];
			a  = new double[dim];
			
			x0  = new double[dim];
			xp0 = new double[dim];
			
			//copy data
			for (int i=0; i<dim; i++) {
				 m[i] =  mm[i];
				sd[i] = ssd[i];
				 a[i] =  aa[i];
				
				//The two initial points need to be in the support
				//ANY point ("within the same galaxy of the main bulk of the objective")
				// commonly works.
				x0[i]  = a[i] + (m[i] - a[i]) * Un01(); //Unab( a[i], m[i]); //Uniform dist.
				xp0[i] = a[i] + (m[i] - a[i]) * Un01(); //Unab( a[i], m[i]);
				
				//It is a very good idea to have a random way to generate initial points.
				//ANY silly distribution may be used.
				//Indeed, this allows comparing different initial points. 
					
			}
			
			
		}
		
		//We need to provide two initial points to the twalk
		double *Getx0() { return x0; }
		double *Getxp0() { return xp0; }
		
		//This is our destructor
		~TruncIndNor() {
			
			delete m;
			delete sd;
			delete a;
			
			delete x0;
			delete xp0;
		}
		
		
		//This method is requiered by obj_fcn
		void show_descrip() const { printf("Obejctive function: Truncated independent Normals.\n"); }


		//This is also requiered by obj_fnc.  A pionter to an array of doubles is passed and
		//if these are within the support of the objective insupport returns 1, otherwise returns 0.
		int insupport(double *x) {
		
			//We check the support:
			//NB: I use gsl's fcmp to compare double's.  See the prototype in ranfun.h .
			for (int i=0; i<get_dim(); i++) //get_dim() is inhiereted from obj_fcn 
				if (fcmp( a[i], x[i]) == 1) //x[i] < a[i] 
					return 0;
			
			return 1;
		
		}
		
		//This is the last method required, the actual objective function.
		//it is guarateed that right after calling insupport and only if it returns 1,
		//eval is called.
		//Then some pretreatment may be done in insupport.  For example, if
		//a correlation parameter is used, the resulting correlation matrix should be
		//positive define.  The support can only be done once positive definess is checked
		//but the cholesky factorization or what ever treatment can be then used in eval
		//without needing to to redo everything.
		
		//It is then assumed then that x is within the support
		//prime=0 means it is x and prime=1 means it is x' (see the twalk paper)
		//eval return -log of the objective.
		double eval(double *x, int prime) {
		
			double U=0;
			for (int i=0; i<get_dim(); i++) 
				U += 0.5*sqr( (x[i] - m[i])/sd[i] ); //sqr is defined in ranfun.h
			
			return U;
		}
		
		
		//We inhierit AccPars(int prime) and NotAcc(int prime) that are called right after eval
		//and tells you what happend with x.  They are defined as virtual so they can be overloaded
		//to do some postprocessing or cleaning.  We ignore them here. 
		
		
		
			
		//************ NOW SEE example.c  *******************

};



#endif
