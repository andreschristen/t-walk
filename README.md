# t-walk
Various implementations of the t-walk: C, C++, R, Python, Matlab and Julia

The t-walk is a "A General Purpose Sampling Algorithm for Continuous Distributions" to sample from many objective functions (specially suited for posterior distributions using non-standard models that would make the use of common algorithms and software difficult); it is an MCMC that does not required tuning.  However, as mentioned in the paper, it may not perform well in some examples and fine tuned samplers to specific objective densities should perform better than the t-walk.

It is now implemented in Python, R, C++, C (native stand alone), MatLab and Julia; see the directories.

The paper is:

Christen, J.A. and Fox, C. (2010), "A General Purpose Sampling Algorithm for
Continuous Distributions (the t-walk)", Bayesian Analysis, 5(2), 263-282.
http://ba.stat.cmu.edu/journal/2010/vol05/issue02/christen.pdf

The t-walk is routinely used running quasi-unsupervised in some aplications (like rbacon, an R package) with more than 1,000 users.  

