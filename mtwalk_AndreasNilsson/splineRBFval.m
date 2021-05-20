function y = splineRBFval(c,alpha,ab,x)
% evaluate spline RBF at points in vector x
% RBF is defined by centers c, coeffs alpha and poly ab (ax + b)

nc = length(c);
nx = length(x);

alphaN = NSplineProj(c,alpha); % project to natural spline

y = ((abs(c(:,ones(1,nx))' - x(:,ones(1,nc)))).^3)*alphaN + polyval(ab,x);
 