function [x,res,it,time]=Steepest_Descent(A,b,x0,tol,itmax)
% STEEPEST DESCENT Resolution of linear systems
%   [x,res]=Steepest_Descent(A,b,x0,tol,itmax) solves Ax = b by
%   Steepest Descent using x0 as initial guess. The algorithm
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
    tol=sqrt(eps);
    itmax=10000;
end

x = x0;
r = b - A*x;
res(1) = norm(b-A*x,2);
it = 1; time = cputime();
while(it <= itmax && res(it)>tol)
    alpha = (r'*r)/(r'*A*r);
    x = x + alpha*r;
    r = b - A*x;
    res(it+1)=norm(b-A*x,2);
    it=it+1;
end
time = cputime()-time;
end

