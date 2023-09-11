function [x,res,iter] = CG(A,b,x0,tol,maxiter)

% CG conjugate gradient method
% [X,res,iter] = CG(A,b,x0,n) computes an approximation for the solution
% of the linear system Ax = b performing n steps of the conjugate
% gradient method starting with x0, and returns in res the norm of the residual

% normb = norm( b );
% if  ( normb == 0.0 )
%     normb = 1.0;
% end
% residual
r = b - A*x0;
x = x0;
% res(1) = norm(r)/normb;
res(1) = norm(r);
p = zeros(size(b));
rho = 1;
it = 1;
while(it<=maxiter && res(it)>tol)
    rhoold = rho; rho = r'*r;
    beta = rho/rhoold;
    p = r + beta*p;
    z = A*p;
    alpha = rho/(p'*z);
    x = x + alpha*p;
    r = r - alpha*z;
%     res(it+1) = norm(r)/normb;
    res(it+1) = norm(r);
    %pause
    it=it+1;
end
iter = it;






