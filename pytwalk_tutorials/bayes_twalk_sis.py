#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 19 10:00:27 2021

@author: jac

Bayesina inference using a SIS model for epidemic data.
The example uses covid19 data from Mexico city.

Let $y_i$ be dayly reports of COVID19, representing observations of $I(t_i)-I(t_{i-1})$
where $\frac{dI}{dt} = \beta I (N - I)$ with $I(0) = I_0$, $t_1 <  t_2 < \cdots < t_n$.
A Poisson model is proposed such as
$$
y_i = Po(I(t_i)-I(t_{i-1})), i=1,2,\ldots,n.
$$
The SIS ODE model above has the following analytical solution
$$
I(t) = \frac{N I_0}{(N - I_0)e^{-\beta t} + I_0}.
$$
$N$ is the city size, and is take as known.
The interest is to estimate $\beta$ and $I_0$, with priors
$$
\beta \sim Gamma( a_0, b_0) ~~\text{and}~~ I_0 \sim Gamma( a_1, b_1, loc), I_0 < N.
$$

"""

import os

from pickle import dump, load

from numpy import array, exp, log, diff, append, arange, mean, quantile, zeros
from scipy.stats import gamma 
from matplotlib.pyplot import subplots

from pytwalk import pytwalk



class sis:
    def __init__( self, name, N, data, times, times1, times_pred, times_pred1):
        """Save data and a name (string of characters) for the model."""
        self.name = name #model instance name
        self.y = data
        self.times = times
        self.times1 = times1
        self. times_pred = times_pred
        self. times_pred1 = times_pred1
        self.N = N
        self.npars = 2 #example ... this is the number of parameters in the model
        self.par_names = [r'\beta', r'I_0'] #example
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

    def PlotData(self, title=""):
        """Plots the data."""
        fig, ax = subplots()
        ax.plot( self.times, self.y, '-')
        ax.plot( self.times, self.y, 'o')
        ax.set_title(title)
        ax.set_xlabel("day")
        ax.set_ylabel("Registered cases")
        fig.tight_layout()

    def I( self, theta, t, N):
        """The analytical solution for I(t), vectorized."""
        b, I0 = theta
        return (N*I0)/((N - I0)*exp(-b*t) + I0)

    def loglikelihood( self, theta):
        ### Evaluate I(t), t=times, the Poisson intensities are the differences. 
        la = diff( self.I( theta, t=self.times1, N=self.N)) + 0.001 #Trick to avoid la=0
        return sum(-la + self.y*log(la))

    def logprior( self, theta, a0=2, b0=0.1, a1=1.0, b1=10.0, loc=1.0):
        beta , I0 = theta   
        return gamma.logpdf( beta, a0, scale=b0) + gamma.logpdf( I0, a1, scale=b1, loc=loc)

    def Supp( self, theta):
        """Support function:
           returns True if the parameters theta are in the support, False otherwise.
        """
        beta , I0 = theta        
        rt = (0 < beta)
        rt *= (1 < I0) and (I0 < N)
        return rt

    def Energy( self, x ):
        """-log of the posterior, input for the twalk."""
        theta = x # Translate x into parameters
        return -(self.loglikelihood( theta )  + self.logprior(theta))

    def Inittheta( self, a0=2, b0=0.1, a1=1.0, b1=10.0, loc=1.0):
        """Simulates initial values for the twalk, simulates from the prior."""
        return array([gamma.rvs( a0, scale=b0), gamma.rvs( a1, scale=b1, loc=loc)])
    
    def RunMCMC( self, T=500000, burnin=1000):
        """Runs the twalk with T iterations and does an output analysis with burnin(=1000)."""
        self.twalk.Run( T=T, x0=self.Inittheta(), xp0=self.Inittheta())
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

    def PlotMCMC(self):
        """Costum plots here.  Example plot of the first 4 parameters."""
        fig, ax = subplots( nrows=2, ncols=2)
        
        ### Plot the burnin removed and thinned with the iat t-walk log-post trace.
        ax[0,0].plot( -self.samples[:,-1], '-')
        ax[0,0].set_xlabel(r"Iteration")
        ax[0,0].set_ylabel(r"Log-post")

        ### Plot marginal posterior for beta
        ax[0,1].hist( self.samples[:,0], density=True)
        ax[0,1].set_xlabel(r"$\beta$")

        ### Plot marginal posterior for I_0
        ax[1,0].hist( self.samples[:,1], density=True)
        ax[1,0].set_xlabel(r"$I_0$")
        
        ### Plot the fit and the predictive
        beta_pm, I0_pm, dummy = mean( self.samples, axis=0)
        ax[1,1].plot( self.times_pred, diff(self.I( [ beta_pm, I0_pm], t=self.times_pred1, N=self.N )), 'r-')
        solns = zeros((self.essize,len(self.times_pred)))
        for i,theta in enumerate(self.samples[:,0:2]):
            solns[i,:] = diff(self.I( self.samples[i,0:2], t=self.times_pred1, N=self.N ))
        qunt = zeros(( len(self.times_pred), 5))
        for i in range(len(self.times_pred)):
            qunt[i,:] = quantile( solns[:,i], q=[0.1, 0.25, 0.5, 0.75, 0.9])
        ax[1,1].fill_between( self.times_pred, qunt[:, 0], qunt[:, 4], color='blue', alpha=0.25)
        ax[1,1].fill_between( self.times_pred, qunt[:, 1], qunt[:, 3], color='blue', alpha=0.25)
        ax[1,1].plot( self.times_pred, qunt[:, 2], 'r-')                

        ax[1,1].plot( self.times, y, '-')
        ax[1,1].plot( self.times, y, '.')
        ax[1,1].set_xlabel("day")
        ax[1,1].set_ylabel("Registered and predicted cases")
        
        fig.tight_layout()
        

        
if __name__ == '__main__':
    
    ### COVID19 cummulative reportes cases, starting 2 March 2020
    covid_vmx_datos = array([ 2., 3., 3., 3., 3., 3., 5., 6., 6., 6., 9., 10., 15., 23., 32., 33., 36., 43., 47.,
                   61., 80., 84., 90., 110., 131., 179., 264., 317., 343., 385., 455., 500., 563., 651.,
                   766., 872., 1045., 1211., 1310., 1512., 1617., 1853., 1997.])

    y = append( covid_vmx_datos[0], diff(covid_vmx_datos)) #dayly cases
    N = 20e6     # Population in the valley of Mexico (Mexico city metro area)
    times = arange(0,len(y))    #dayly cases 
    times1 = arange(0,len(y)+1) 
    times_pred = arange(0,len(y)+14) # Two week prediction 
    times_pred1 = arange(0,len(y)+15) # Two week prediction 

    covid_vmx = sis( name="covid_vmx", N=20e6, data=y, times=times, times1=times1, times_pred=times_pred, times_pred1=times_pred1)
    covid_vmx.PlotData(title="COVID19 dayly cases, Vally of Mexico\nRegister from 2 March 2020")
    #covid_vmx.RunMCMC()
    #covid_vmx.PlotMCMC()
    
    
        
