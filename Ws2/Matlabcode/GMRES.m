%************************************************
function [x,error,iter,time] = GMRES(A,b,x,m,itmax,epsi)
%************************************************

% gmres.m solves the linear system Ax=b
% using the Generalized Minimal residual ( GMRESm ) method with restarts .
%
% input   A        REAL nonsymmetric positive definite matrix
%         x        REAL initial guess vector
%         b        REAL right hand side vector
%         m        INTEGER number of iterations between restarts
%         itmax   INTEGER maximum number of iterations
%         epsi      REAL error tolerance
%
% output  x        REAL solution vector
%         error    REAL error norm
%         iter     INTEGER number of iterations performed

% initialization

time = cputime();
normb = norm(b);
if  ( normb == 0.0 )
    normb = 1.0;
end
% residual
r = b - A*x;
error(1) = norm(r)/normb;
if ( error(1) < epsi )
    return;
end
[n,n] = size(A);
V(1:n,1:m+1) = zeros(n,m+1);   % Vm+1=[Vm|qm+1]
H(1:m+1,1:m) = zeros(m+1,m);   % Hessenberg matrix
cs(1:m) = zeros(m,1);
sn(1:m) = zeros(m,1);

e1    = zeros(n,1);            % basic vector
e1(1) = 1.0;
iter = 1;  % step of iterator
itt = 1;
while iter <= itmax            % begin iteration
    r = b-A*x;
    V(:,1) = r/norm(r);
    s = norm(r)*e1;
    for j = 1:m                % construct orthonormal
        itt = itt+1 ;          % basis using Gram-Schmidt
        w = A*V(:,j);
        for i = 1:j
            H(i,j)= w'*V(:,i);
            w = w - H(i,j)*V(:,i);
        end
        %     size(H)
        H(j+1,j) = norm(w);
        V(:,j+1) = w/H(j+1,j);
        % We tranform the Hessenberg matrix H into a triangular matrix by applying Givens rotation
        for k = 1:j-1                              % apply Givens rotation
            temp     =  cs(k)*H(k,j) + sn(k)*H(k+1,j);
            H(k+1,j) = -sn(k)*H(k,j) + cs(k)*H(k+1,j);
            H(k,j)   = temp;
        end
        cs(j) = H(j,j)/sqrt(H(j,j)^2 + H(j+1,j)^2);
        sn(j) = H(j+1,j)/sqrt(H(j,j)^2 + H(j+1,j)^2);
        temp   = cs(j)*s(j);                        % approximate residual norm
        s(j+1) = -sn(j)*s(j);
        s(j)   = temp;
        H(j,j) = cs(j)*H(j,j) + sn(j)*H(j+1,j);
        H(j+1,j) = 0.0;
        % Which also provides the error at iteration j+1
        
        error(itt+1) = abs(s(j+1)) / normb;
        
        if ( error(itt+1) <= epsi )                    % update approximation
            y = H(1:j,1:j)\s(1:j);                 % and exit
            x = x + V(:,1:j)*y;
            % error(i+1) = abs(s(i+1)) / bnrm2;
            break;
        end
    end
    
    if ( error(itt+1) <= epsi ), break, end
    y = H(1:m,1:m)\s(1:m);
%        x = x + V*y;                                   % update approximation
    x = x + V(:,1:m)*y;                                 % update approximation
    % compute residual
    r = b - A*x;
    s(j+1) = norm(r);
    error(itt+1) = s(j+1) / normb;                        % check convergence
    if ( error(itt+1) <= epsi )
        break;
    end
    iter = iter+1;
end

iter = itt;
time = cputime()-time;