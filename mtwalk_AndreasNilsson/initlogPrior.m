function rmat = initlogPr(c)
% set up Hessian for log prior

nc = length(c);

D2 = 6*abs(c(:,ones(1,nc))' - c(:,ones(1,nc)));

R = zeros(nc,nc);
for count = 1:nc-1
    R(count:count+1,count:count+1) = R(count:count+1,count:count+1) + ...
        [1 1/2;1/2 1]*(c(count+1)-c(count))/3;
end

rmat = D2*R*D2;