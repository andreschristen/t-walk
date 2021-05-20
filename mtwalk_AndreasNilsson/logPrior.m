function logPr = logPrior(alpha,ab,rho,rmat,c)
% evaluate log prior

logPr = -rho*alpha'*rmat*alpha -(rho/10)*((sum(alpha))^2 + (sum(alpha.*c))^2);