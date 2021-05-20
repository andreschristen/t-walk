function [xvec,dvec,sigma] = MockData
% function file to make mock data

sigma = 0.1;
nx = 30;
xvec = linspace(0,1,nx)';
dvec = func1(xvec) + randn(size(xvec))*sigma;

%-------------------------
function y = func(x)
%true function is piecewise linear
y = 1-2*abs(x-1/2);

%-------------------------
function y = func1(x)
%true function is sine bump
y = sin(x*pi);