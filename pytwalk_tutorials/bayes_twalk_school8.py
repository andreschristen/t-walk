#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 19 10:00:27 2021

@author: jac

Example of the "8 scholls", commonly used as the first example in BUGG, STAN etc.
see: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

The proble is this.  There are 8 schools for which a treatment is taken (couching)
and resultson the effect, with their sd deviations, are given.
That is the data below in the school8_data dictionary.

The modelo explained in Section 5.5 of Gelman et al (2003) is this:
    
Let $y_j \pm \sigma_j$ be the treatment effect of school $j$, with its sd $\sigma_j$.
\begin{align*}
y_j | \theta_j, \sigma_j & \sim N( \mu + \tau \eta_j, \sigma_j) ~~i=1,2,\ldots,J\\
\eta_j &\sim N(0,1) .
\end{align*}
Costant improper priors are given for the population mean $\mu$ and the sacle $\tau$.
Although in STAN $\log(\tau)$ is said to be constant (not the same).
Here $\tau$ is constant with its support $\tau >0$.  

"""

import os

from pickle import dump, load

from numpy import array, append
from scipy.stats import norm, uniform
from matplotlib.pyplot import subplots

from pytwalk import pytwalk

### There are other nicer corner, or pairs, plots, use a better one.
from PlotCorner import PlotCorner


### 8 School data 
school8_data = { 'names': ["A","B","C","D","E","F","G","H"],\
            'effect': [28.39,7.94,-2.75,6.82,-.64,.63,18.01,12.16],\
            'sde': [14.9, 10.2, 16.3, 11.0, 9.4, 11.4, 10.4, 17.6]}


### The class.  Start from a copy of the templete, and fill in the blanks.
class school8:
    def __init__( self, name, data):
        self.name = name ### Needs a name for saving the output
        ### Transalte the data into the class data memmbers
        self.y = array(data['effect'])
        self.sd = array(data['sde'])
        self.s_names = array(data['names'])
        self.J = len(self.y) # Number of schools
        ### Total number of parameters
        self.npars = self.J + 2
        ### and their names.
        ### If an array holds the model pameters (e.g. theta, x), they are arrange
        ### always in this orther:
        self.par_names = [r'\mu', r'\tau'] + [r'\eta_%d' % (i) for i in range(self.J)]
        ### Define the t-walk, for the MCMC
        self.twalk = pytwalk( n=self.npars, U=self.Energy, Supp=self.Supp)

        ### Check if file self.name + '_samples.pkl', to read previously generated
        ### MCMC samples.
        if os.path.isfile(self.name + '_samples.pkl'): # samples file exists
            print("File %s with mcmc samples exists, loading samples ..." % (self.name + '_samples.pkl',), end=' ')
            self.samples = load(open(self.name + '_samples.pkl', 'rb'))
            self.essize = self.samples.shape[0]
            print(" done.")
        else:
            print("File %s with mcmc samples does not exist, run RunMCMC first." % (self.name + '_samples.pkl',))
        
        
    def loglikelihood( self, theta):
        ### Translate the parameters:
        mu, tau = theta[0:2]
        eta = theta[2:]
        ### The log-likelihood is the sum of log of the gaussians
        return sum( norm.logpdf( self.y, loc=mu + tau*eta, scale=self.sd))

    def logprior( self, theta):
        mu, tau = theta[0:2]
        eta = theta[2:]
        ### Only the priors for eta_j, mu and tau have constant priors. 
        return sum(norm.logpdf( eta, loc=0, scale=1))        

    def Supp( self, theta):
        """Support function:
           returns True if the parameters theta are in the support, False otherwise.
        """
        mu, tau = theta[0:2]
        eta = theta[2:]
        rt = tau > 0 ### Only tau > 0, the other have no restrictions
        return rt

    def Energy( self, x ):
        """-log of the posterior, input for the twalk."""
        return -(self.loglikelihood( x )  + self.logprior(x))

    def Inittheta( self):
        """Simulates initial values for the twalk, normally from the prior."""
        eta = norm.rvs(size=self.J)
        mu = uniform.rvs(loc=-40, scale=80)
        tau = uniform.rvs(loc=1, scale=20)
        tmp = array([mu, tau])
        return append( tmp, eta)

    def RunMCMC( self, T=500000, burnin=1000):
        """Runs the twalk with T iterations and does an output analysis with burnin(=1000)."""
        self.twalk.Run( T=T, x0=self.Inittheta(), xp0=self.Inittheta())
        self.AnaSaveMCMC( burnin=burnin)
        
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
        """Some custom plots for this inference problem."""
        par_labels = ["$%s$" % (n) for n in self.par_names]
        ### Corner plot of 4 parameters
        PlotCorner( s8.samples[:,0:4], names=par_labels )
        ### Posteriors for mu and tau
        fig, ax = subplots( nrows=1, ncols=2)
        ax[0].hist(self.samples[:,0])
        ax[0].set_xlabel(par_labels[0])
        ax[1].hist(self.samples[:,1])
        ax[1].set_xlabel(par_labels[1])
        fig.tight_layout()
        ### Posterior of all eta_j's, plotted horizontally using box-plots         
        fig, ax = subplots()
        ax.boxplot( self.samples[:,2:-1], vert=False)
        ax.set_yticklabels(self.s_names)
        ax.set_ylabel("Schools")
        ax.set_xlabel("Normalized shifts (posteriors for $\eta_j$'s, plotted with boxplots)")
        fig.tight_layout()
        
if __name__ == '__main__':
    
    s8 = school8( data=school8_data, name="school8")
    s8.RunMCMC()
    s8.PlotMCMC()
    
    
    
    