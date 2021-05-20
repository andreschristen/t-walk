#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 19 10:00:27 2021

@author: jac

Templete for doing Bayesian inference using the t-walk:
"""

import os

from pickle import dump, load

from pytwalk import pytwalk

from PlotCorner import PlotCorner


class model:
    def __init__( self, name, data):
        """Save data and a name (string of characters) for the model."""
        self.name = name #model instance name
        self.y = data
        self.npars = 4 #example ... this is the number of parameters in the model
        self.par_names = [r'\alpha', r'\beta_1', r'\beta_2', 'E'] #example
        ### define the twalk object
        self.twalk = pytwalk( n=self.npars, U=self.Energy, Supp=self.Supp)
        ### Check if self.name + '_samples.pkl' exists, to load previous (thinned) MCMC samples
        if os.path.isfile(self.name + '_samples.pkl'): # samples file exists
            print("File %s with mcmc samples exists, loading samples ..." % (self.name + '_samples.pkl',), end=' ')
            self.samples = load(open(self.name + '_samples.pkl', 'rb'))
            self.essize = self.samples.shape[0]
            print(" done.")
        else:
            print("File %s with mcmc samples does not exist, run RunMCMC first." % (self.name + '_samples.pkl',))
        ###### Need to define the loglikelihood and logprior
        
    def loglikelihood( self, theta):
        pass

    def logprior( self, theta):
        pass

    def Supp( self, theta):
        """Support function:
           returns True if the parameters theta are in the support, False otherwise.
        """
        pass

    def Energy( self, x ):
        """-log of the posterior, input for the twalk."""
        theta = x # Translate x into parameters
        return -(self.loglikelihood( theta )  + self.logprior(theta))

    def Inittheta( self):
        """Simulates initial values for the twalk, normally from the prior."""
        pass
    
    def RunMCMC( self, T=500000, burnin=1000):
        """Runs the twalk with T iterations and does an output analysis with burnin(=1000)."""
        self.pytwalk.Run( T=T, x0=self.Inittheta(), xp0=self.Inittheta())
        self.AnaSaveMCMC(burnin=burnin)

    def AnaSaveMCMC( self, burnin):
        """Output analysis of the MCMC output, with burnin."""
        ### This calculates the Intehrated Autocorrelation Time
        ### provides some more information on the screen
        ### and makes a trace plot of the logpost, after burnin.
        self.iat = int(self.twalk.Ana(start=burnin)[0,0])
        self.burnin = burnin
        ### Uses iat for thining the MCMC and saveing the samples.
        self.samples = self.twalk.Output[burnin::(self.iat),:] ## Burn in and thining
        self.essize = self.samples.shape[0]
        print("\nEffective sample size: %d" % (self.essize,))
        dump( self.samples, open(self.name + '_samples.pkl', 'wb'))

    def PlotMCMC( self, burnin=10000):
        """Costum plots here.  Example plot of the first 4 parameters."""
        par_labels = ["$%s$" % (n) for n in self.par_names]
        ### Corner plot of 4 parameters
        PlotCorner( s8.samples[:,0:4], names=par_labels )
        ### Example of marginal posterior with 4 parameters
        fig, ax = subplots( nrows=2, ncols=2)
        ax[0,0].hist(self.samples[:,0])
        ax[0,0].set_xlabel(par_labels[0])
        ax[0,1].hist(self.samples[:,1])
        ax[0,1].set_xlabel(par_labels[1])
        ax[1,0].hist(self.samples[:,2])
        ax[1,0].set_xlabel(par_labels[2])
        ax[1,1].hist(self.samples[:,3])
        ax[1,1].set_xlabel(par_labels[3])
        

        
if __name__ == '__main__':
    
    mo = model( data=, name="...")
    mo.RunMCMC(T=100000, burnin=1000)
    
    
        
