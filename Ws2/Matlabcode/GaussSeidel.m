function [x,res,it,time]=GaussSeidel(A,b,x0,tol,itmax)
% GAUSSSEIDEL Resolution of linear systems by Gauss-Seidel
%   [x,res]=GAUSSSEIDEL(A,b,x0,tol,itmax) solves Ax = b by
%   Gauss-Seidel using x0 as initial guess. The algorithm
%   stops when maxiter operations have been performed, or when the L2 norm of
%   residual becomes smaller than tol. The output x is the approximation of the
%   solution, and res is a vector containing the history of the residuals.

if(nargin < 3)
    error('Not enough arguments!!');
end
if(nargin < 5)
    itmax=10000;
end
if(nargin <4)
    tol=sqrT(eps);
    itmax=10000;
end

D = diag(diag(A));
E = -tril(A) + D;
F = -triu(A) + D;

x = x0;
res(1) = norm(b-A*x,2);
it = 1; time = cputime();
while(it <= itmax && res(it)>tol)
    %while(i<=maxiter) % to see the behavior of the iterations
    x = inv(D-E)*(b + F*x);
    res(it+1)=norm(b-A*x,2);
    it=it+1;
end
time = cputime()-time;
end

