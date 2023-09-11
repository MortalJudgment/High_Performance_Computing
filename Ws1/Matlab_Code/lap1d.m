function A=lap1d(n)
% lap1d one dimensional finite difference approximation
%   A=lap1d(n) computes a sparse finite difference
%   approximation of the one dimensional operator -Delta on the
%   domain Omega=(0,1) using n interior points

h = 1/(n+1);
e = ones(n,1);
A=spdiags([-e/h^2 2/h^2*e -e/h^2],[-1 0 -1],n,n);
