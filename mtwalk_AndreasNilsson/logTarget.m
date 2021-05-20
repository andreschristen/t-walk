function lf = logTarget(x,p)
% calculate log target distribution for use in t-walk

%[rho rmat c xvec dvec sigma] = param;

alpha = x(1:end-2);
ab = x(end-1:end);

logP = logPrior(alpha,ab,p.rho,p.rmat,p.c);
logL = logLikelihood(p.c,alpha,ab,p.xvec,p.dvec,p.sigma);

lf = logP + logL;

