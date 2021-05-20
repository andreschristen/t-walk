function alphaN = NSplineProj(c,alpha)
% project coefficients onto space of natural splines

nc = length(c);

gamma = -inv([nc sum(c); sum(c) sum(c.^2)])*[sum(alpha); sum(alpha.*c)];

alphaN = (alpha  + gamma(2)*c) + gamma(1);
