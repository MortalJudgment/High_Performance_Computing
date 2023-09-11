function [x,res,it,time] = CG(A,b,x0,tol,itmax)
% CG conjugate gradient method
% [X,R,a,b]=CG(A,b,x0,n) computes an approximation for the solution
% of the linear system Ax=b performing n steps of the conjugate
% gradient method starting with x0, and returns in res the norm of the residual
x = x0; r = b-A*x;
res(1) = norm(r,2);
p = zeros(size(b));
rho = 1;
it = 1; time = cputime();
while(it<=itmax && res(it)>tol)
    rhoold = rho; rho = r'*r;
    beta = rho/rhoold;
    p = r + beta*p;
    z = A*p;
    alpha = rho/(p'*z);
    x = x + alpha*p;
    r = r - alpha*z;
    res(it+1) = norm(r,2);
    it = it+1;
end
time = cputime()-time;

