function [x,res,i]=CG(A,b,x0,maxiter,tol);

% CG conjugate gradient method
% [X,R,a,b]=CG(A,b,x0,n) computes an approximation for the solution
% of the linear system Ax=b performing n steps of the conjugate
% gradient method starting with x0, and returns in res the norm of the residual
x=x0; r=b-A*x;
res(1)=norm(r,inf);
p=zeros(size(b));
rho=1;
i=1;

normb = norm( b );
if  ( normb == 0.0 ),
    normb = 1.0;
end
% residual
r=b-A*x;
res(1) = norm( r )/normb;

while(i<=maxiter & res(i)>tol)
    rhoold=rho; rho=r'*r;
    beta=rho/rhoold;
    p=r+beta*p;
    z=A*p;
    alpha=rho/(p'*z);
    x=x+alpha*p;
    r=r-alpha*z;
    res(i+1)= norm(r) / normb ;
    %pause
    i=i+1;
end
iter=i;






