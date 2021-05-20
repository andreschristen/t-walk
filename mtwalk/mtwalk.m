function [xxp,lt,acc] = twalk_p(logTarget,nsamp,x0,xp0,p)
% Function file implements t-walk sampling from target distribution
%
% [xxp,lt,acc] = twalkparam(logTarget,nsamp,x0,xp0,p)
%
% Input arguments:
% logTarget name of function that calculates the log of the target density.
%           Should return -Inf for points not in support.
% nsamp     number of steps
% x0,xp0    starting values (n by 1 column vectors, not equal)
% p         parameter to pass to logtarget function
% Output arguments:
% xxp       2*n by (nsamp+1) array of output samples
%           xxp(1:n,:) is chain of x values
%           xxp(n+(1:n),:) is chain of xp values
% lt        chain of values of log target density
% acc       vector of acceptance rates for each move, and mean value.
             
% Colin Fox, abril 2006
% modified 17 May 2009 to pass a parameter to the logtarget function

% Andreas Nilsson, July 2016
% changed parameters aw and at according to Christen and Fox (2010)
% introduced [phi] to only move proportion of the parameters

% Initial values
x = x0;
xp = xp0;

ltx = feval(logTarget,x,p);
ltxp = feval(logTarget,xp,p);

if (ltx == -Inf) || (ltxp == -Inf)
    disp('Initial values out of support:  '), eval('x'), eval('xp')
    return
end

MoveRatio = [0 60 60 1 1]; % first ratio for move 0 -- i.e. don't move.
MoveProb = cumsum(MoveRatio(1:end-1)/sum(MoveRatio));

n1 = 4; % expected number of parameters to be moved
pphi = min(length(x),n1)/length(x);

prop = zeros(size(MoveRatio));       % collect number of proposals
acc = zeros(size(MoveRatio));        % collect acceptance ratios
xxp = zeros(2*length(x),nsamp+1);    % to store output
lt = zeros(2,nsamp+1);               % store values of log target density

xxp(:,1) = [x; xp];
lt(:,1) = [ltx;ltxp];

for iter = 1:nsamp
    dir = rand < 0.5;                    % choose a direction
    kernum = sum(MoveProb <= rand);      % choose a kernel (0 to #Moves-1)
    phi = rand(length(x),1)<pphi;        % set parameters to move
    prop(kernum+1) = prop(kernum+1) + 1;

    switch kernum
        case 0 % do nothing, ensures aperiodicity ... not needed (see paper)
            if dir, yp = xp; else y = x; end, lgr = 0; % always accepted
        case 1 % traverse move
            if dir, [yp,lgr] = Simh1(xp,x,phi); else [y,lgr] = Simh1(x,xp,phi); end
        case 2 % walk move
            if dir, [yp,lgr] = Simh2(xp,x,phi); else [y,lgr] = Simh2(x,xp,phi); end
        case 3 % hop move
            if dir, [yp,lgr] = Simh3(xp,x,phi); else [y,lgr] = Simh3(x,xp,phi); end
        case 4 % blow move
            if dir, [yp,lgr] = Simh4(xp,x,phi); else [y,lgr] = Simh4(x,xp,phi); end
    end
    
    if dir
        y = x; lty = ltx;
        ltyp = feval(logTarget,yp,p);
        lalpha = ltyp - ltxp + lgr;
    else
        yp = xp; ltyp = ltxp;
        lty = feval(logTarget,y,p);
        lalpha = lty - ltx + lgr;
    end
        
    if (rand < exp(lalpha))
        % accepted
        acc(kernum+1) =  acc(kernum+1) + sum(phi)/length(x);
        
        x = y; ltx = lty;
        xp = yp; ltxp = ltyp;
    end
    xxp(:,iter+1) = [x;xp];
    lt(:,iter+1) = [ltx;ltxp];
end

% calculate acceptance rates
acc = [acc./(prop+eps), sum(acc)/nsamp];

%-------------------------------------------------------------------------
function beta = psi1
% function beta = psi1
% return beta distributed as psi1

at = 6; % PARAMETER_at in distribution

if rand < (at-1)/(2*at)
    beta = rand^(1/(1+at));
else
    beta = rand^(1/(1-at));
end

%-------------------------------------------------------------------------
function [h,lgr] = Simh1(x,xp,phi)
% [h,lgr] = Simh1(x,xp)
% simulate traverse kernel

beta = psi1; % simulate psi1 to give beta
h = x;
h(phi) = xp(phi) + beta*(xp(phi) - x(phi));

nI = sum(phi);
lgr = (nI-2)*log(beta); % since psi1(beta) = psi1(1/beta)

%-------------------------------------------------------------------------
function [h,lgr] = Simh2(x,xp,phi)
% [h,lgr] = Simh2(x,xp)
% walk kernel

aw = 1.5; % PARAMETER_aw in distribution

u = rand(size(x));
z = (aw/(1+aw))*(-1 + 2*u + aw*u.^2); % z ~ 1/sqrt(1+z) in [-a/(1+a),a]

h =  x + phi.*(x - xp).*z;

lgr = 0; % As z has density 1/sqrt(1+z)

%-------------------------------------------------------------------------
function [h,lgr] = Simh3(x,xp,phi)
% [h,lgr] = Simh3(x,xp)
% hop kernel

if sum(phi)>0
    sigma = max(abs(xp(phi) - x(phi)));
else
    sigma = 0;
end
h = x + phi.*sigma/3.*randn(size(x)); % iid standard normal

% evaluate log of ratio of densities g(x|h,x')/g(h|x,x')
lgr = lg3(x,h,xp,phi) - lg3(h,x,xp,phi);

%-------------------------------------------------------------------------
function lg = lg3(h,x,xp,phi)
% log density for hop move, up to constant not affecting ratio

if sum(phi)>0
    sigma = max(abs(xp(phi) - x(phi)));
else
    sigma = 0;
end
nI = sum(phi);
lg = -nI*log(sigma) - sum((h-x).^2)/(2*sigma^2);

%-------------------------------------------------------------------------
function [h,lgr] = Simh4(x,xp,phi)
% [h,lgr] = Simh4(x,xp)
% blow kernel

if sum(phi)>0
    sigma = max(abs(xp(phi) - x(phi)));
else
    sigma = 0;
end
h = x;
h(phi) = xp(phi) + sigma*randn(size(x(phi))); % iid standard normal

% evaluate log of ratio of densities g(x|h,x')/g(h|x,x')
lgr = lg4(x,h,xp,phi) - lg4(h,x,xp,phi);

%-------------------------------------------------------------------------
function lg = lg4(h,x,xp,phi)
% log density for blow move, up to constant not affecting ratio

if sum(phi)>0
    sigma = max(abs(xp(phi) - x(phi)));
else
    sigma = 0;
end
nI = sum(phi);
lg = -nI*log(sigma) - sum((h-xp).^2)/(2*sigma^2);

