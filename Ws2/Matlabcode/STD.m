function[x,error,iter]=STD(A,b,x,itmax,epsi)  

%************************************************

% STD.m solves the linear system Ax=b with the Steepest Descent method
%
% input   A        REAL nonsymmetric positive definite matrix
%         x        REAL initial guess vector
%         b        REAL right hand side vector
%         itmax   INTEGER maximum number of iterations
%         epsi      REAL error tolerance
%
% output  x        REAL solution vector
%         error    REAL error norm
%         iter     INTEGER number of iterations performed



% initialization
normb = norm( b );
if  ( normb == 0.0 ),
    normb = 1.0;
end
% residual
r=b-A*x;
error(1) = norm( r )/normb;

if ( error(1) < epsi ) 
    return; 
end

iter=1;  % step of iterator
while iter <= itmax  
    % begin iteration
    Ar = A*r ;
    rr = dot(r,r) ;
    rar = dot(r,Ar) ;
    alp = rr/rar ;
    x = x + alp*r  ;   
    % compute residual
    r=r -alp*Ar;
    error(iter+1) = norm(r) / normb   ;                   % check convergence
     if ( error(iter+1) <= epsi )
        break;
    end
    iter=iter+1;
end

