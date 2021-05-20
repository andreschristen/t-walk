#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#  twalktutorial.py
#  
#  Examples for the twalk implementation in Python,
#  Created by J Andres Christen, jac at cimat.mx .
#
#
"""


from numpy import ones, zeros, log, array
from numpy.random import uniform
from matplotlib.pyplot import figure, close

from pytwalk import pytwalk

close('all')

"""
########################################################
### This sets an MCMC for n independent normals ########
### Three basic ingredients are needed:
### 1) define -log of the objective (posterior)
### 2) define the support function, return True if parameters are in the support
### 3) Define *two* initial points for the t-walk.
### Each parameter needs to be different in the initial points.
### Commonly, we include a function to simulae two random points, within the support is used.

### Commonly (not always) two initial points remotely within the same 'galaxy'
### as the actual 'effective' support
### (where nearly all of the probability is) of the objectuve are only needed.
"""

def NormU(x):
    """Defines the 'Energy', -log of the objective, besides an arbitrary constant."""
    return sum(0.5*x**2)

def NormSupp(x):
    return True

Nor = pytwalk( n=10, U=NormU, Supp=NormSupp)

### This runs the twalk, with initial points 10*ones(10) and zeros(10)
### with T=100000 iterations
print("\n\nRunning the dimension 10 Gaussian example:""")
Nor.Run( T=100000, x0=10*ones(10), xp0=zeros(10))

### This does a basic output analysis, with burnin=1000, in particular the IAT:
figure(1)
Nor.Ana(start=1000)
figure(2)
Nor.Hist( par=0, start=1000)

"""
########################################################
######### Example of a product of exponentilas #########
### Objective function:
### $f(x_1,x_2,x_3,x_4,x_5)=\prod_{i=1}^5 \lambda exp(-x_i \lambda)$
"""

lambdas = [ 1., 2., 3., 4., 5.]

def ExpU(x):
	"""-log of a product of exponentials"""
	return sum(x * lambdas)

def ExpSupp(x):
	return all(0 < x)

Exp = pytwalk( n=5, U=ExpU, Supp=ExpSupp)

### This runs the twalk, with initial points ones(40) and 40*ones(40)
### with T=50000 iterations
print("\n\nRunning the product of exponentials example:""")
Exp.Run( T=50000, x0=30*ones(5), xp0=40*ones(5))
figure(3)
Exp.Ana(start=2000)
figure(4)
Nor.Hist( par=2, start=1000) # Plot the third parameter


"""
#########################################################################
##### A more complex example:
##### Related Bernoulli trials #####
##### Suppose x_{i,j} ~ Be( theta_j ), i=0,1,2,...,n_j-1, ind. j=0,1,2
##### But it is known that 0 <  theta_0 < theta_3 < theta_2 < 1 
"""

theta = array([ 0.4, 0.5, 0.7 ])  ### True thetas
n = array([ 20, 15, 40]) ### sample sizes
#### Simulated synthetic data: 
r = zeros(3)
for j in range(3):
    ### This su¡imulates each Bernoulli: uniform(size=n[j]) < theta[j]
    ### but we only need the sum of 1's (sum of True's)
	r[j] = sum(uniform(size=n[j]) < theta[j])

### Defines the support.  This is basically the prior, uniform in this support
def ReBeSupp(theta):
	rt = True
	rt &= (0 < theta[0])
	rt &= (theta[0] < theta[1])
	rt &= (theta[1] < theta[2])
	rt &= (theta[2] < 1)
	return rt

#### It is a good idea to have a function that produces random initial points,
#### indeed always withinh the support
#### In this case we simulate from something similar to the prior
def ReBeInit():
	theta = zeros(3)
	theta[0] = uniform( low=0, high=1)
	theta[1] = uniform( low=theta[0], high=1)
	theta[2] = uniform( low=theta[1], high=1)
	return theta

####### The function U (Energy): -log of the posterior
def ReBeU(theta):
	return -1*sum(r * log(theta) + (n-r)*log(1.0-theta))
	
###### Define the twalk instance
ReBe = pytwalk( n=3, U=ReBeU, Supp=ReBeSupp)

#### This runs the twalk
print("\n\nRunning the related Bernoullis example:""")
ReBe.Run( T=100000, x0=ReBeInit(), xp0=ReBeInit())

### First plot a trace of the log-post, to identify the burn-in
figure(5)
ReBe.TS() 

### This will do a basic output analysis, with burnin=5000
figure(6)
ReBe.Ana(start=5000)

### And some histograms
figure(7)
ReBe.Hist( par=0 )
figure(8)
ReBe.Hist( par=1 )
figure(9)
ReBe.Hist( par=2 )	

### Then save it to a text file, with each column for each paramter
### plus the U's in the last column, that is T+1 rows and n+1 colums.
### This may in turn be loaded by other programs
### for more sophisticated output analysis (eg. BOA).
ReBe.Save("RelatedBer.txt")

### You may access the (T+1) X (n+1) output matrix directly with
#ReBe.Output

###### All methods should have help lines: ReBe.Hist?















