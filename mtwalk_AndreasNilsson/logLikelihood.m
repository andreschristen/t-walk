function logL = logLikelihood(c,alpha,ab,xvec,dvec,sigma)
% evaluate log Likelihood

logL = -sum((dvec - splineRBFval(c,alpha,ab,xvec)).^2)/(2*sigma^2);