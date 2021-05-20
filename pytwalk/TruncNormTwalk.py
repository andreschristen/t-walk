#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 11:25:43 2019

@author: jac
"""

from numpy import array, var, mean, sqrt
from scipy.stats import norm

from pytwalk import pytwalk

n = 5
#m  = array([2.0, 3.9, 5.2, 1.0, 7.8]) ##means
#sd = array([1.2, 2.1, 3.3, 0.1, 2.3]) ##std. dev's
#a  = array([0.0, 0.0, 2.0, 0.0, 3.0]) ##Lower boundaries

m  = array([2.0]*5) ##means
sd = array([1.0]*5) ##std. dev's
a  = array([-100.00]*5) ##Lower boundaries



def TruncNormSupp(x):
    return all(a < x);

def TruncNormU(x):
    return sum(0.5*((x - m)/sd)**2)

def TruncNormInit(n=n):
    x = a.copy()
    while not(TruncNormSupp(x)):
        x = m + sd*norm.rvs(size=n)
    return x

TruncNorm = pytwalk( n=n, U=TruncNormU, Supp=TruncNormSupp)
                    #ww=[0.0000, 0.4918, 0.4918, 0.0082+0.0082, 0.0])

TruncNorm.Run( T=1000000, x0=TruncNormInit(), xp0=TruncNormInit())

print("Mean:\n", [mean(TruncNorm.Output[:,i]) for i in range(5)] )

print("Sd dev:\n", [sqrt(var(TruncNorm.Output[:,i])) for i in range(5)] )

"""
the Hop move now works ... Felipe found the bugg in line 531

WHY THE HOP MOVE DOES NOT WORK ...

WHY IT was not used in cpptwalk ??


"""


