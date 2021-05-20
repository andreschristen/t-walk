% script file to set up spline RBF and sample posteior using twalk

nc = 11;
cmin = 0; cmax = 1;
c = linspace(cmin,cmax,nc)'; % linearly spaced centres

[xvec,dvec,sigma] = MockData;

alpha = zeros(size(c)); ab = [0;0]; % zero initial spline

rmat = initlogPrior(c);
rho = 1;
logP = logPrior(alpha,ab,rho,rmat,c);
logL = logLikelihood(c,alpha,ab,xvec,dvec,sigma);

x0 = zeros(nc+2,1); % zero starting state
xp0 = x0+1;

p = struct('rho', rho, 'rmat', rmat, 'c', c, 'xvec', xvec, 'dvec', dvec, 'sigma', sigma);

hold off
for count = 1:1000
%     [xxp,ltxp,acc] = twalkparam('logTarget',100,x0,xp0,p);
%     [xxp,ltxp,acc] = twalk_c('logTarget',100,x0,xp0,p);
    [xxp,ltxp,acc] = twalk_p('logTarget',100,x0,xp0,p);
    alphahat=xxp(1:nc,end); abhat=xxp(nc+1:nc+2,end);
    if rem(count,10)==0
    plot(xvec,dvec,'ro',xvec,splineRBFval(c,alphahat,abhat,xvec),'b',c,splineRBFval(c,alphahat,abhat,c),'bo')
    axis([0 1 -.2 1.2])
    drawnow
    end
    x0 = xxp(1:nc+2,end); 
    xp0 = xxp(nc+3:end,end);
end