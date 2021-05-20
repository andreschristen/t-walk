/*
   This code is derived from the C++ version of t-Walk "cpptwalk-beta-1.0"
   Feb 2008 version, kindly supplied by J. Andres Christen.

   Apart from converting this to a stand-alone C program I can claim no
   other credit.  I freely admit that C++ is a better language for this
   code because of the ability to derive the application dependent
   Objective Function class from an abstract base class.

   The t-Walk code is an implementation of some original and clever ideas
   for a universal Metropolis-Hastings MCMC sampler in the paper "A General
   Purpose Sampling Algorithm for Continuous Distributions (the t-walk)" by
   J. Andres Christen and Colin Fox, published in Bayesian Analysis (2010)
   Vol 5, Number 2, pp. 263-282.  The paper is well worth a read!

   Compile with cc Ctwalk.c -o Ctwalk -lm

   Tony Begg <tony.begg@dataventures.com>
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define New(objects, have, TYPE)                                            \
   (objects) = (TYPE *)malloc((have) * sizeof(*(objects)));                 \
   if ((objects) == NULL)                                                   \
   {                                                                        \
      fprintf(stderr, "New(): Out of Heap Memory (%d)\n", have);            \
      abort();                                                              \
      exit(1);                                                              \
   }

/*
   fcmp
   Copyright (c) 1998-2000 Theodore C. Belding
   University of Michigan Center for the Study of Complex Systems
   <mailto:Ted.Belding@umich.edu>
   <http://www-personal.umich.edu/~streak/>
  
   This file is part of the fcmp distribution. fcmp is free software;
   you can redistribute and modify it under the terms of the GNU Library
   General Public License (LGPL), version 2 or later.
*/
int fcmp(double x1, double x2)
{
   static const double epsilon = 1.0e-11; /* used for t-walk */
   int exponent;
   double delta;
   double difference;
   
   /*
      Get exponent(max(fabs(x1), fabs(x2))) and store it in exponent.
   */
   frexp(fabs(x1) > fabs(x2) ? x1 : x2, &exponent);
   /*
      Do the comparison.
      delta = epsilon * pow(2, exponent)
   */
   delta = ldexp(epsilon, exponent); 
   difference = x1 - x2;
   if (difference > delta)
   {
      return 1; /* x1 > x2 */
   }
   else if (difference < -delta) 
   {
      return -1;  /* x1 < x2 */
   }
   else /* -delta <= difference <= delta */
   {
      return 0;  /* x1 == x2 */
   } 
}

/*
   Random is based on rng/taus.c from the GNU Scientific Library.
   Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
   
   NormalDev is based on randist/gauss.c from the GNU Scientific Library.
   Copyright (C) 1996, 1997, 1998, 1999, 2000, 2006 James Theiler, Brian Gough
   Copyright (C) 2006 Charles Karney

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or (at
   your option) any later version.
   
   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
struct Random
{
   unsigned long int s1, s2, s3;
};

typedef struct Random Random;

static unsigned long Random32(Random *rng)
{
#define MASK 0xffffffffUL
#define TAUSWORTHE(s,a,b,c,d) (((s & c) << d) & MASK) ^ ((((s << a) & MASK)^s) >> b)

   rng->s1 = TAUSWORTHE (rng->s1, 13, 19, 4294967294UL, 12);
   rng->s2 = TAUSWORTHE (rng->s2, 2, 25, 4294967288UL, 4);
   rng->s3 = TAUSWORTHE (rng->s3, 3, 11, 4294967280UL, 17);
   return (rng->s1 ^ rng->s2 ^ rng->s3);
}

/*
   Uniform real random number between 0 and 1.
*/
double Un01(Random *rng)
{
   return (double)Random32(rng) / 4294967296.0;
}

/*
   Uniform real random number between values a and b where b > a.
*/
double Unab(Random *rng, double a, double b)
{
   return (((b - a) * ((double)Random32(rng) / 4294967296.0)) + a);
}

void RandomSeed(Random *rng, unsigned long int s)
{
   if (s == 0) s = 1;      /* default seed is 1 */
#define LCG(n) ((69069 * n) & 0xffffffffUL)
   rng->s1 = LCG(s);
   if (rng->s1 < 2) rng->s1 += 2UL;
   rng->s2 = LCG(rng->s1);
   if (rng->s2 < 8) rng->s2 += 8UL;
   rng->s3 = LCG(rng->s2);
   if (rng->s3 < 16) rng->s3 += 16UL;
   /* "warm it up" */
   Random32(rng);
   Random32(rng);
   Random32(rng);
   Random32(rng);
   Random32(rng);
   Random32(rng);
   return;
}

double NormalDev(Random *rng, double mean, double sigma)
{
   double u, v, x, y, Q;
   const double s = 0.449871;  /* Constants from Leva */
   const double t = -0.386595;
   const double a = 0.19600;
   const double b = 0.25472;
   const double r1 = 0.27597;
   const double r2 = 0.27846;

   do /* This loop is executed 1.369 times on average  */
   {
      /*
         Generate a point P = (u, v) uniform in a rectangle enclosing
         the K+M region v^2 <= - 4 u^2 log(u).
         u in (0, 1] to avoid singularity at u = 0.
      */
      u = 1.0 - Un01(rng);
      /*
         v is in the asymmetric interval [-0.5, 0.5).  However v = -0.5
         is rejected in the last part of the while clause.  The
         resulting normal deviate is strictly symmetric about 0
         (provided that v is symmetric once v = -0.5 is excluded).
      */
      v = Un01(rng) - 0.5;
      /*
         Constant 1.7156 > sqrt(8/e) (for accuracy); but not by too
         much (for efficiency).
      */
      v *= 1.7156;
      /*
         Compute Leva's quadratic form Q.
      */
      x = u - s;
      y = fabs(v) - t;
      Q = x * x + y * (a * y - b * x);
      /*
         Accept P if Q < r1 (Leva)
         Reject P if Q > r2 (Leva)
         Accept if v^2 <= -4 u^2 log(u) (K+M)
         This final test is executed 0.012 times on average.
      */
   }
   while (Q >= r1 && (Q > r2 || v * v > -4 * u * u * log (u)));
   return mean + (sigma * (v / u)); /* Return slope plus mean */
}

/*
   Here is our random number generator object.  It is used for
   all our deviate needs.
*/
static Random RNG;

/*
   These are the parameter settings a_w, a_t and n1, as defined in the paper.
*/
#define PARAMETER_aw (1.5)
#define PARAMETER_at (6.0)
/*
   Expected number of coordinates to be moved, parameter n_1.
*/
#define EXP_MOV_COOR (4)

#define min(x,y) ((x) < (y) ? (x) : (y))

/*
   Some useful functions for vectors of doubles.
*/
/*
   Allocate storage for a vector of length n.
*/  
double *Vector(int n)
{
   double *v;

   New(v, n, double);
   return v;
}

/*
   Free storage for a vector.
*/  
void FreeVector(double *v)
{
   free((void *)v);
}

/*
   Copy vector src to vector dest.
*/
void VectorCopy(double *src, double *dest, int n)
{
   int i;

   for (i = 0; i < n; i++) dest[i] = src[i];
}

/*
   Subtract vector v2 from v1, elementwise, and store result in diff.
*/
void VectorDiff(double *v1, double *v2, int n, double *diff)
{
   int i;

   for (i = 0; i < n; i++) diff[i] = v1[i] - v2[i];
}

/*
   Find the maximum absolute value of a vector (masked by phi) and
   return its index.
*/
void VectorIndexMax(double *v, int n, int *ix, int *phi)
{
   int i, max_i;

   max_i = 0;
   for (i = 0; i < n; i++)
   {
      max_i = ((fcmp((double)phi[max_i] * fabs(v[max_i]),
                     (double)phi[i] * fabs(v[i])) == -1) ? i : max_i);
   }
   *ix = max_i;
}

/*
   Compare vector v1 with v2, elementwise, and return 1 if all
   elements agree (are equal), else return 0.
*/
int VectorCmp(double *v1, double *v2, int n)
{
   int i = 0;

   while ((fcmp(v1[i], v2[i]) == 0) && (i < n)) i++;
   if (i == n) return 1;
   return 0;
}

/*
   Print out a vector on the file pointer fp in a particular format.
   Note that large n will lead to very long lines.
*/ 
void VectorPrint(FILE *fp, double *v, int n)
{
   int i;

   fprintf(fp, "\n");
   for (i = 0; i < n; i++)
   {
      fprintf(fp, "\t%13.6g", v[i]);
   }
}

/*
   The ObjFunc object is used to define the Posterior Probability being
   explored by the t-walk MCMC.  The methods are therefore problem
   specific and user defined.  Example methods from the C++ code example
   are used below for illustration.  Here in this comment we give the
   basic interface (in C++ the virtual base class from which user ObjFunc
   objects are derived). 

   C   struct ObjFunc
   C   {
   C      int dim;
   C   };
   C   
   C   typedef struct ObjFunc ObjFunc;
   C   
   C   void ObjFuncCTOR(ObjFunc *S, int dim)
   C   {
   C      S->dim = dim;
   C   }
   C   
   C   void ObjFuncDTOR(ObjFunc *S)
   C   {
   C   }
   C   
   C   int ObjFuncGetDim(ObjFunc *S)
   C   {
   C      return S->dim;
   C   }
   C   
   C   void ObjFuncAccPars(ObjFunc *S, int prime)
   C   {
   C   }
   C   
   C   void ObjFuncNotAcc(ObjFunc *S, int prime)
   C   {
   C   }
   C   
   C   void ObjFuncShowDescrip(ObjFunc *S)
   C   {
   C   }
   C   
   C   int ObjFuncInSupport(ObjFunc *S, double *x)
   C   {
   C   }
   C   
   C   double ObjFuncEval(ObjFunc *S, double *x, int prime)
   C   {
   C   }
*/
/*
   There now follows J. Andres Christen's TruncIndNor example objective
   function object distributed with the twalk package, adapted here
   for the C language.
*/
struct ObjFunc
{
   int dim; /* number of parameters in model */
   double *m; /* pointer to the set of mean's for each normal */
   double *sd; /* pointer to the array of standard dev's for each normal */
   double *a; /* the lower boundaries for the support of each normal */
   double *x0, *xp0; /* two initial points */
};

typedef struct ObjFunc ObjFunc;

/*
   Our constructor.  The ObjFunc constructor takes one integer, dim, as
   argument, which is the dimension of the objective function.  We take
   that as well and pass it on.
*/
void ObjFuncCTOR(ObjFunc *S, int dim, double *m, double *sd, double *a)
{
   int i;

   S->dim = dim;
   /*
      Allocate space.
   */
   New(S->m, dim, double);
   New(S->sd, dim, double);
   New(S->a, dim, double);
   New(S->x0, dim, double);
   New(S->xp0, dim, double);
   /*
      Copy data from the caller.
   */
   for (i = 0; i < dim; i++)
   {
      S->m[i] =  m[i];
      S->sd[i] = sd[i];
      S->a[i] =  a[i];
      /*
         The two initial points need to be in the support.
         ANY point "within the same galaxy of the main bulk of the
         objective" commonly works.
      */
      S->x0[i] = Unab(&RNG, a[i], m[i]);
      S->xp0[i] = Unab(&RNG, a[i], m[i]);
   }
}

/*
   This is our destructor.
*/
void ObjFuncDTOR(ObjFunc *S)
{
   free((void *)S->m);
   free((void *)S->sd);
   free((void *)S->a);
   free((void *)S->x0);
   free((void *)S->xp0);
}

/*
   We need to provide two initial points to the twalk.
   Note these are extra methods to the base class ObjFunc but
   are called to provide constructor arguments to the TWalk
   object, so do not appear within the common generic code.
*/
double *ObjFuncGetx0(ObjFunc *S)
{
   return S->x0;
}

double *ObjFuncGetxp0(ObjFunc *S)
{
   return S->xp0;
}

int ObjFuncGetDim(ObjFunc *S)
{
   return S->dim;
}

/*
   We inherit AccPars(int prime) and NotAcc(int prime) that are called
   right after Eval and tell us what happened with x.  They are defined
   as virtual so they can be overloaded to do some postprocessing or
   cleaning.  We ignore them here. 
*/
void ObjFuncAccPars(ObjFunc *S, int prime)
{
}

void ObjFuncNotAcc(ObjFunc *S, int prime)
{
}

/*
   This is required by ObjFunc.
*/
void ObjFuncShowDescrip(ObjFunc *S)
{
   fprintf(stdout, "Objective function: Truncated Independent Normals.\n");
}

/*
   This is also required by ObjFunc.
   A pointer to an array of doubles is passed and if these are within
   the support of the objective function InSupport returns 1, otherwise
   returns 0.
*/
int ObjFuncInSupport(ObjFunc *S, double *x)
{
   int i;

   /*
      We check the support:
      NB: We use Theodore C. Belding's fcmp (also in GNU Scientific Library)
      to compare doubles.
   */
   for (i = 0; i < S->dim; i++)
   {
      /*
         Test x[i] < a[i].
      */
      if (fcmp(S->a[i], x[i]) == 1) return 0; 
   }
   return 1;
}

/*
   This is the last method required, the actual objective function.
   It is guaranteed that Eval is called right after calling InSupport
   but only if it returns 1.
   Therefore, some pretreatment may be done in InSupport.  For example,
   if a correlation parameter is used, the resulting correlation matrix
   should be positive definite.  The support can only be done once
   positive definiteness is checked but the cholesky factorization or
   whatever treatment can be then used in Eval without needing to redo
   everything.
      
   It is then assumed that x is within the support prime=0 means it is
   x and prime=1 means it is x' (see the twalk paper).

   Eval returns minus the log of the objective.
*/
double ObjFuncEval(ObjFunc *S, double *x, int prime)
{
   int i;
   double tmp;
   double U = 0.0;

   for (i = 0; i < S->dim; i++) 
   {
      tmp = (x[i] - S->m[i]) / S->sd[i];
      U += 0.5 * tmp * tmp;
   }
   return U;
}
/*
   End of example Objective Function code.
   The rest of the code, except main() is generic.
*/

/*
   -log f_a(beta).
   It is assumed that beta > 0
*/
double Ufbeta(double beta, double a)
{
   double b = 0.0;

   if (beta < 1.0)
   {
      b = -a * log(beta);
   }
   if (beta >= 1.0)
   {
      b = a * log(beta);
   }
   return (-log(a - 1.0) - log(a + 1.0) + log(2.0 * a) + b);
}

double Simfbeta(double a)
{
   if (Un01(&RNG) < (a - 1.0) / (2.0 * a))
   {
      return (exp(1.0 / (a + 1.0) * log(Un01(&RNG))));
   }
   else
   {
      return (exp(1.0 / (1.0 - a) * log(Un01(&RNG))));
   }
}

/*
   The argument is parameter aw in the paper.
*/
double Phi2Sim(double aw)
{
   double u = Un01(&RNG);

   return ((aw / (1.0 + aw)) * ((2.0 * u) + (aw * u * u) - 1.0));
}

/*
   Kernel object.  We do not use overloading of methods so need a kernel
   type identifier.  I do not think K0 is used elsewhere in the code.
*/
#define K0 (0)
#define K1 (1)
#define K2 (2)
#define K3 (3)
#define K4 (4)

struct Kernel
{
   int type; /* to distinguish between K0, K1, K2, K3 and K4 */
   double *h; /* auxiliary vector */
   double sigma;
   double *rest;
};

typedef struct Kernel Kernel;

void KernelCTOR(Kernel *S, int type)
{
   S->type = type;
   S->h = NULL;
}

void KernelDTOR(Kernel *S)
{
}

void KernelSetH(Kernel *S, double *h)
{
   S->h = h;
}

double *KernelSimH(Kernel *S, double *x, double *xp,
                   int n, double beta, int *phi)
{
   int i, j;
   double *retn;

   switch (S->type)
   {
      case K0:
         retn = x;
         break;

      case K1:
         for (i = 0; i < n; i++)
         {
            if (phi[i] == 1)
            {
               S->h[i] = xp[i] + (beta * (xp[i] - x[i]));
            }
            else
            {
               S->h[i] = x[i];
            }
         }
         retn = S->h;
         break;

      case K2:
         for(i = 0; i < n; i++)
         {
            S->h[i] =
               x[i] + ((double)phi[i] * (x[i] - xp[i]) * Phi2Sim(PARAMETER_aw));
         }
         retn = S->h;
         break;

      case K3:
         VectorDiff(xp, x, n, S->rest);
         VectorIndexMax(S->rest, n, &i, phi);
         S->sigma = fabs(S->rest[i]) / 3.0;
         for (j = 0; j < n; j++)
         {
            S->h[j] =
               x[j] + ((double)phi[j] * S->sigma * NormalDev(&RNG, 0.0, 1.0));
         }
         retn = S->h;
         break;

      case K4:
         VectorDiff(xp, x, n, S->rest);
         VectorIndexMax(S->rest, n, &i, phi);
         S->sigma = fabs(S->rest[i]);
         for (j = 0; j < n; j++)
         {
            if (phi[j] == 1)
            {
               S->h[j] = xp[j] + (S->sigma * NormalDev(&RNG, 0.0, 1.0));
            }
            else
            {
               S->h[j] = x[j];
            }
         }
         retn = S->h;
         break;

      default:
         retn = NULL;
         break;
   }
   return retn;
}

double KernelGU(Kernel *S, double *h, double *x, double *xp, int n)
{
   int i, j;
   double intProd, retn;

   switch (S->type)
   {
      case K0:
         if (VectorCmp(h, x, n))
         {
            retn = 1.0;
         }
         else
         {
            retn = 0.0;
         }
         break;

      case K1:
         retn = -2.0;
         break;

      case K2:
         retn = 1.0;
         break;

      case K3:
         if (! VectorCmp(x, xp, n))
         {
            intProd = 0.0;
            for (j = 0; j < n; j++)
            {
               intProd += (h[j] - x[j]) * (h[j] - x[j]);
            }
            /*
               It is assumed that SimH has just been called and
               we have the correct sigma.
            */
            retn = ((((double)n * 0.5) * log(2.0 * M_PI)) +
                    ((double)n * log(S->sigma)) +
                    (0.5 * (1.0 / (S->sigma * S->sigma)) * intProd));
         }      
         else
         {
            retn = -1.0;
         }
         break;

      case K4:
         intProd = 0.0;
         for(j = 0; j < n; j++)
         {
            intProd += (h[j] - xp[j]) * (h[j] - xp[j]);
         }
         /*
            It is assumed that SimH has just been called and
            we have the correct sigma.
         */
         retn = ((((double)n * 0.5) * log(2.0 * M_PI)) +
                 ((double)n * log(S->sigma)) +
                 (0.5 * (1.0 / (S->sigma * S->sigma)) * intProd));
         break;

      default:
         retn = 0.0;
         break;
   }
   return retn;
}

double Kernelfbeta(Kernel *S, double beta, double a)
{
   double b = 0.0;

   if (0.0 <= beta && beta < 1.0)
   {
      b = pow(beta, a);
   }
   if (beta > 1.0)
   {
      b = pow(1.0 / beta, a);
   }
   return (((a - 1.0) * (a + 1.0) * b) / (2.0 * a));
}

double KernelFbeta(Kernel *S, double beta, double a)
{
   double b1 = 0.0;
   double b2 = 0.0;

   if (0.0 <= beta && beta <= 1.0)
   {
      b1 = ((a - 1.0) / (2.0 * a)) * pow(beta, a + 1.0);
   }
   if (beta > 1.0)
   {
      b1 = (a - 1.0) / (2.0 * a);
      b2 = ((a + 1.0) / (2.0 * a)) * (1.0 - pow(1.0 / beta, a - 1.0));
   }
   return (b1 + b2);
}

void KernelSetRest(Kernel *S, double *rest)
{
   S->rest = rest;
}

/*
   The t-walk object.
*/
struct TWalk
{
   /*
      Pointer to objective function object.
   */
   ObjFunc *Obj;
   /*
      Initial values and the pool of points updated in each iteration
   */
   double *x;
   double *xp;
   double U, Up;
   double *h;
   double *rest;
   /*
      Dimension of the objective function domain
   */
   int n; /* x, xp dimension */
   /*
      Parameters needed by methods Simulation and OneMove
   */
   double acc; /* acceptance counter */
   int val;
   Kernel *ker;
   double propU, propUp, *y, *yp, dir,  W1, W2, A, aux, beta;
   int *phi;
   double mapU, *mapx;
   /*
      Saving scheme: if save_every < 0 only accepted iterations.
      if > 0 only it %% abs(save_every) == 0 (1 => all iterations).
      if 0, save all iterations and save the file recacc with kernel
      acceptance rates etc., for debugging purposes.
   */
   int save_every;
   int debugg;
   /*
      Transition kernels
   */
   Kernel k1;
   Kernel k2;
   Kernel k3;
   Kernel k4;
   int nphi;
   double pphi;
};

typedef struct TWalk TWalk;

/*
   Selection between (y,h(x,xp)) and (h(xp,x),yp)
*/
double TWalkSelectPivot(TWalk *S)
{
   double aux = Un01(&RNG);

   return aux;
}

/*
   Selection between (y,h(x,xp)) and (h(xp,x),yp)
*/
int *TWalkSelectPhi(TWalk *S, int *phi)
{
   int i;

   S->nphi = 0;
   for (i = 0; i < S->n; i++) 
   {
      if (Un01(&RNG) < S->pphi)
      {
         phi[i] = 1;
         S->nphi++;
      }
      else
      {
         phi[i] = 0;
      }
   }
   return phi;
}

void TWalkCTOR(TWalk *S, ObjFunc *Obj, double *x, double *xp, int n)
{
   S->Obj = Obj;
   S->x = x;
   S->xp = xp;
   S->n = n;
   KernelCTOR(&S->k1, K1);
   KernelCTOR(&S->k2, K2);
   KernelCTOR(&S->k3, K3);
   KernelCTOR(&S->k4, K4);
   S->h = Vector(n);
   S->rest = Vector(n);
   KernelSetH(&S->k1, S->h);
   KernelSetH(&S->k2, S->h);
   KernelSetH(&S->k3, S->h); KernelSetRest(&S->k3, S->rest);
   KernelSetH(&S->k4, S->h); KernelSetRest(&S->k4, S->rest);
   S->pphi = (double)min(n, EXP_MOV_COOR) / (double)n;
   S->mapx = Vector(n);
   S->y = Vector(n);
   S->yp = Vector(n);
   S->acc = 0.0;
   S->A = 0.0;
   New(S->phi, n, int);
}

void TWalkDTOR(TWalk *S)
{
   FreeVector(S->h);
   FreeVector(S->rest);
   FreeVector(S->mapx);
   FreeVector(S->y);
   FreeVector(S->yp);
   free((void *)S->phi);
}

Kernel *TWalkSelectKernel(TWalk *S, int *type)
{
   double aux = Un01(&RNG);

   if (0.0000 <= aux && aux < 0.0082) /* Kernel K3, hop */
   {
      *type = 3;
      return &S->k3;
   }
   else if (0.0082 <= aux && aux < 0.0164) /* kernel K4, blow */
   {
      *type = 4;
      return &S->k4;
   }
   else if (0.0164 <= aux && aux < 0.5082) /* kernel K1, traverse */
   {
      *type = 1;
      return &S->k1;
   }
   else if (0.5082 <= aux && aux < 1.0) /* kernel K2, walk */
   {
      *type = 2;
      return &S->k2;
   }
   *type = 2; /* should not reach here */
   return &S->k2;
}

/*
   Initialization.
*/
int TWalkInit(TWalk *S, double *x, double *xp)
{
   int i;
   int n = S->n;

   if (x != NULL)
   {
      VectorCopy(x, S->x, n);
   }
   if (xp != NULL)
   {
      VectorCopy(xp, S->xp, n);
   }
   if (ObjFuncInSupport(S->Obj, S->x))
   {
      S->U = ObjFuncEval(S->Obj, S->x, 0);
      ObjFuncAccPars(S->Obj, 0);
   }
   else
   {
      fprintf(stdout, "twalk: parameters x out of support:");
      for (i = 0; i < n; i++)
      {
         fprintf(stdout, "\t%11.6g\n", S->x[i]);
      }
      fprintf(stdout, "\n");
      return 0;
   }
   if (ObjFuncInSupport(S->Obj, S->xp))
   {
      S->Up = ObjFuncEval(S->Obj, S->xp, 1);
      ObjFuncAccPars(S->Obj, 1);
   }
   else
   {
      fprintf(stdout, "twalk: parameters xp out of support:");
      for (i = 0; i < n; i++)
      {
         fprintf(stdout, "\t%11.6g\n", S->xp[i]);
      }
      fprintf(stdout, "\n");
      return 0;
   }
   S->mapU = S->U;
   VectorCopy(S->x, S->mapx, n);
   S->propU = S->U;
   S->propUp = S->Up;
   VectorCopy(S->x, S->y, n);
   VectorCopy(S->xp, S->yp, n);
   return 1;
}

int TWalkOneMove(TWalk *S)
{
   int n = S->n;

   /*
      Selection of the pivot and transition kernel
   */   
   S->ker = TWalkSelectKernel(S, &S->val);
   S->dir = TWalkSelectPivot(S);
   (void)TWalkSelectPhi(S, S->phi);
   if (S->dir < 0.5)
   {
      /*
         x is the pivot.
         Beta is a dummy parameter for some kernels.
      */   
      S->beta = Simfbeta(PARAMETER_at);
      /*
         yp is the proposal
      */
      VectorCopy(KernelSimH(S->ker, S->xp, S->x, n, S->beta, S->phi), S->yp, n);
      VectorCopy(S->x, S->y, n);
      S->propU = S->U;
      /*
         Verify that the proposal is in the obj function domain
      */
      if (ObjFuncInSupport(S->Obj, S->yp))
      {
         /*
            Evaluate the obj function in the proposal.
         */   
         S->propUp = ObjFuncEval(S->Obj, S->yp, 1);
         /*
            Compute the acceptance probability.
         */
         S->W1 = KernelGU(S->ker, S->yp, S->xp, S->x, n);
         S->W2 = KernelGU(S->ker, S->xp, S->yp, S->x, n);
         if (S->W1 == -1.0 || S->W2 == -1.0)
         {
            S->A = 0.0;
         }
         else if (S->W1 == -2.0 && S->W2 == -2.0)
         {
            S->A = exp((S->U - S->propU) + (S->Up - S->propUp) +
                       ((double)(S->nphi - 2) * log(S->beta)));
         }
         else
         {
            S->A = exp((S->U - S->propU) + (S->Up - S->propUp) + (S->W1 - S->W2));
         }
      }
      else
      {
         S->A = 0.0;
      }
   }
   else
   {
      /*
         xp is the pivot.
         Repeating the procedure above but using xp as pivot.
      */
      S->beta = Simfbeta(PARAMETER_at);
      VectorCopy(KernelSimH(S->ker, S->x, S->xp, n, S->beta, S->phi), S->y, n);
      VectorCopy(S->xp, S->yp, n);
      S->propUp = S->Up;
      if (ObjFuncInSupport(S->Obj, S->y))
      {
         S->propU = ObjFuncEval(S->Obj, S->y, 0);
         S->W1 = KernelGU(S->ker, S->y, S->x, S->xp, n);
         S->W2 = KernelGU(S->ker, S->x, S->y, S->xp, n);
         if (S->W1 == -1.0 || S->W2 == -1.0)
         {
            S->A = 0.0;
         }
         else if (S->W1 == -2.0 && S->W2 == -2.0)
         {
            S->A = exp((S->U - S->propU) + (S->Up - S->propUp) +
                       ((double)(S->nphi - 2) * log(S->beta)));
         }
         else
         {
            S->A = exp((S->U - S->propU) + (S->Up - S->propUp) + (S->W1 - S->W2));
         }
      }
      else
      {
         S->A = 0.0;
      }
   }
   S->aux = Un01(&RNG);
   if (S->aux < S->A)
   {
      /*
         Accepted
      */
      S->acc += (double)S->nphi / (double)n; /* proportion of moved parameters */
      VectorCopy(S->y, S->x, n);
      S->U = S->propU;
      VectorCopy(S->yp, S->xp, n);
      S->Up = S->propUp;
      if (S->dir >= 0.5)
      {
         /*
            y is accepted.
         */
         ObjFuncAccPars(S->Obj, 0);
         if (fcmp(S->U, S->mapU) == -1)
         {
            /*
               U < mapU
            */
            S->mapU = S->U;
            VectorCopy(S->x, S->mapx, n);
         }
         return 1; /* y accepted */
      }
      else
      {
         /*
            yp is accepted
         */
         ObjFuncAccPars(S->Obj, 1);
         return -1; /* yp accepted */
      }
   }
   else
   {
      if (S->dir >= 0.5)
      {
         /*
            y is not accepted.
         */
         ObjFuncNotAcc(S->Obj, 0);
      }
      else
      {
         /*
            yp is not accepted.
         */
         ObjFuncNotAcc(S->Obj, 1);
      }
      return 0;
   }
}

/*
   Information messages.
   total iterations, current iteration, start time, current time.
*/
void Remain(int Tr, int it, long sec1, long sec2)
{
   long ax;

   /*
      how many seconds remaining.
   */
   ax = (long)((double)(Tr - it) * ((double)(sec2 - sec1)/(double)it));
   if (ax == 0)
   {
      fprintf(stdout, "\n");
      fflush(stdout);
      return;
   }
   if (ax < 60)
   {
      fprintf(stdout, "Finish in approx. %ld seconds.\n", ax);
      fflush(stdout);
      return;
   }
   if (ax <= 360)
   {
      fprintf(stdout, "Finish in approx. %ld minutes and %ld seconds.\n",
              ax / 60, ax % 60);
      fflush(stdout);
      return;
   }
   if (ax > 360)
   {
      ax += sec2;  /* end time = current time plus seconds remaining */
      fprintf(stdout, "Finish by %s\n", ctime(&ax));
      fflush(stdout);
      return;
   }
}

/*
   Here is the implementation of the control structure part of
   the algorithm.

   We pass in a file pointer for parameter output rather than a
   filename and options to the routine fopen.  This is the only
   functionality change from the original C++ code.  It has the
   disadvantage that the output filename cannot be incorporated
   in informative log messages.

   If TWalkSimulation is called with NULL arguments xx and xxp then
   a simulation can be continued from a previous state, provided
   the TWalk object has not been destroyed in the meantime.
   NULL arguments are also used if x and xp have just been set
   via TWalkCTOR.
*/
int TWalkSimulation(TWalk *S, int iters, FILE *fptr,
                    int save_every, double *xx, double *xxp)
{
   static const int wait = 30;
   FILE *recacc;
   long sec, sec2, ax;
   int n, it, j, j1, rt, acc_it, abs_save_every;

   n = S->n;
   sec = time(NULL); /* begining of the twalk */
   fprintf(stdout, "twalk: %10d iterations to run, %s", iters, ctime(&sec));
   /*
      Initialization.
   */
   if (TWalkInit(S, xx, xxp) == 0)
   {
      exit(1);
   }
   /*
      send an estimation for the duration of the sampling if 
      evaluating the objective function twice takes more than
      one second
   */
   sec2 = time(NULL); /* last time we sent a message */
   fprintf(stdout, "       ");
   Remain(iters, 2, sec, sec2);
   S->save_every = save_every;
   abs_save_every = save_every < 0 ? -save_every: save_every;
   if (S->save_every == 0)
   {
      S->save_every = 1;
      S->debugg = 1;
   }
   else
   {
      S->debugg = 0;
   }
   if (S->debugg)
   {
      if ((recacc = fopen("recacc.dat", "w")) == NULL)
      {
         fprintf(stderr, "Could not open file %s for writing\n", "recacc.dat");
         exit(1);
      }
      fprintf(stdout,
   "twalk: Kernel acceptance rates information to be saved in file recacc.dat\n");
   }
   VectorPrint(fptr, S->x, n);
   fprintf(fptr, "\t%lf", S->U);
   if (S->save_every < 0)
   {
      fprintf(stdout,
        "twalk: Only every %d accepted iterations will be saved in output file\n",
              abs_save_every);
   }
   else
   {
      fprintf(stdout,
         "twalk: Every %d iterations to be saved in output file\n",
              S->save_every);
   }
   j1 = 1; j = 0; acc_it = 0;
   for (it = 1; it <= iters; it++)
   {  
      rt = TWalkOneMove(S);
      if ((rt == 1) || (rt == -1))
      {
         acc_it++;
         if (S->save_every < 0) /* Only accepted iterations are saved */
         {
            if ((acc_it % abs_save_every) == 0)
            {
               VectorPrint(fptr, S->x, n);
               fprintf(fptr, "\t%13.6g", S->U);
            }
         }
         if (S->debugg)
         {
            fprintf(recacc, "%d\t%f\n", S->val, (double)S->nphi / (double)n);
         }
      }
      else
      {
         /*
            Proposal not accepted.
         */
         if (S->debugg)
         {
            fprintf(recacc, "%d\t%f\n", S->val, 0.0);
         }
      }
      if (S->save_every > 0)
      {
         /*
            whether accepted or not, iterations are saved.
         */
         if ((it % S->save_every) == 0)
         {
            VectorPrint(fptr, S->x, n);
            fprintf(fptr, "\t%13.6g", S->U);
         }
      }
      if ((it % (1 << j1)) == 0)
      {
         j1++;
         j1 = min(j1, 10);
         /*
            check the time at least every 2^10 = 1024 iterations
         */
         if (((ax = time(NULL)) - sec2) > ((1 << j) * wait))
         {
            fprintf(stdout, "twalk: %10d iterations so far. ", it);
            Remain(iters, it, sec, ax);
            sec2 = ax;
            j++;
            j1--; /* check the time as often */ 
         }
      }
   }
   if (S->debugg)
   {
      fclose(recacc);
   }
   /*
      return current points
   */
   if (xx != NULL) VectorCopy(S->x, xx, n);
   if (xxp != NULL) VectorCopy(S->xp, xxp, n);
   sec = time(NULL);
   fprintf(stdout,
      "twalk: %4.1f%% of moved parameters per iteration ratio (%.0f / %d)\n",
           100.0*(S->acc / (double)iters), S->acc, iters);
   fprintf(stdout, "twalk: Finished %s\n", ctime(&sec));
   return (int)(S->acc + 0.5);
}

/*
   This is the example main program from the twalk distribution
   for Truncated Independent Normals.
*/
int main(int argc, char **argv)
{
   ObjFunc obj;
   TWalk twalk;
   FILE *fp;
   /*
      Provide all the parameters:
   */
   int dim = 5;
   double m[5] = { 2.0, 3.9, 5.2, 1.0, 7.8 }; /* means */
   double sd[5] = { 1.2, 2.1, 3.3, 0.1, 2.3 }; /* standard deviations */
   double a[5] = { 0.0, 0.0, 2.0, 0.0, 3.0 }; /* Lower (truncation) boundaries */
   
   if (argc != 2)
   {
      fprintf(stderr, "Usage: %s output_file\n", argv[0]);
      exit(1);
   }
   fp = fopen(argv[1], "w");
   if (fp == NULL)
   {
      fprintf(stderr, "Cannot open output file (%s) for writing\n", argv[1]);
      exit(1);
   }
   /*
      This sets the seed for the random number generator, defined as
      a static structure early in this file.  It must be set before
      calling the Truncated Independent Normal ObjFunc constructor
      since it uses Unab uniform random number generator.
   */
   RandomSeed(&RNG, (unsigned long int)53);
   /*
      Open a Truncated Independent Normal object.
   */
   ObjFuncCTOR(&obj, dim, m, sd, a);
   /*
      We now open a twalk object.  Its arguments are a pointer to
      the Objective Function, the initial points x and x' and the
      Objective Function dimension.
   */
   TWalkCTOR(&twalk, &obj, ObjFuncGetx0(&obj),
             ObjFuncGetxp0(&obj), ObjFuncGetDim(&obj));
   /*
      The twalk will use the initial values AND the space provided by x0
      and xp0.  Therefore at the end of the run, we have access to such
      values.
   
      And now we run the twalk. We specify 100000 iterations and a file
      pointer to save the output to.  The 1 is the "save_every" parameter
      which specifies the iterations saving scheme:
         if > 0, save only iterations such that (iters % save_every == 0)
                 (thus 1 => all iterations).
         if < 0, save only accepted iterations such that
                 (iters % abs(save_every) == 0). 
         if = 0, save all iterations, and save debugging information on each
                 kernel's acceptance rates in the file recacc.dat.
   */
   (void)TWalkSimulation(&twalk, 100000, fp, 1, NULL, NULL);
   /*
      The output file from twalk will have (dim + 1) columns (6 in this
      case), one for each parameter, plus a last column with the
      corresponding value of Eval (minus the log of the objective).
      This is very useful for convergence checking.  Only the x and
      not the x' are saved.
   */
   exit(0);
}
