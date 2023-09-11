function[x,error,iter]=SOR(A,b,x,itmax,epsi,omg)  

%************************************************

% GS.m solves the linear system Ax=b with the Gauss Seidel method
%
% input   A        REAL nonsymmetric positive definite matrix
%         x        REAL initial guess vector
%         b        REAL right hand side vector
%         itmax   INTEGER maximum number of iterations
%         epsi      REAL error tolerance
%         omg       REAL  facteur de relaxation
%
% output  x        REAL solution vector
%         error    REAL error norm
%         iter     INTEGER number of iterations performed



% initialization
dd = diag(A) ;
D =  diag( dd ) ;
E = -tril(A) +D ;
F = -triu(A) +D ;
DOE = D - omg*E  ;
Dm1 = inv(DOE);
dm1b = omg*Dm1*b ;
G = Dm1*( omg*F + (1-omg)*D ) ;


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
while iter <= itmax                              % begin iteration
    x = G*x + dm1b ;   
    % compute residual
    r=b-A*x;
    error(iter+1) = norm(r) / normb ;                       % check convergence
     if ( error(iter+1) <= epsi )
        break;
    end
    iter=iter+1;
end

