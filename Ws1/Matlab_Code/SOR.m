function [x,res,it,time]=SOR(A,b,relax,x0,tol,itmax)
% SOR Resolution of linear systems by overrelaxation method
%   [x,res]=SOR(A,b,relax,x0,tol,itmax) solves Ax=b by
%   SOR with parameter relax  using x0 as initial guess. The algorithm
%   stops when maxiter operations have been performed, or when the L2 norm of
%   residual becomes smaller than tol. The output x is the approximation of the
%   solution, and res is a vector containing the history of the residuals.

if(nargin < 4)
    error('Not enough arguments!!');
end
if(nargin < 6)
    itmax = 10000;
end
if(nargin <5)
    tol = sqrT(eps);
    itmax = 10000;
end

E = tril(A,-1);
F = triu(A,1);
D = diag(diag(A));

x = x0;
res(1)=norm(b-A*x,2);
it=1; time = cputime();
while it<=itmax && res(it)>tol
    x = (D/relax+E)\(((1-relax)/relax*D-F)*x+b);
    res(it+1) = norm(b-A*x,2);
    it = it+1;
end
time = cputime()-time;
end

