function [x, error, iter] = Right_PRECGMRES( A,b,y,M,m,itmax,epsi)
% PRECGMRES.m solves the preconditioner linear system
% The left preconditioner solve MAx = b
% The right preconditioner solve AMy = b, then x = My
% Using the above gmres.m, modify it to precgmres.m so as to be able to use
% a right preconditioner M that can be changed at each restart of gmres.
% input   A        REAL nonsymmetric positive definite matrix
%         x        REAL initial guess vector
%         b        REAL right hand side vector
%         M        REAL preconditioner matrix
%         numb_m   INTEGER number of iterations between restarts
%         max_it   INTEGER maximum number of iterations
%         tol      REAL error tolerance
%
% output  x        REAL solution vector
%         error    REAL error norm
%         iter     INTEGER number of iterations performed
% initialization

normb = norm(b);
if  ( normb == 0.0 )
    normb = 1.0;
end

% residual
r = b - A*M*y;
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
iter = 1;                      % step of iterator
itt = 1;
while iter <= itmax            % begin iteration
    r = b - A*M*y;
    V(:,1) = r/norm(r);
    s = norm(r)*e1;
    for j = 1:m                 % construct orthonormal
        itt = itt + 1;          % basis using Gram-Schmidt
        w = A*M*V(:,j);
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
            z = H(1:j,1:j)\s(1:j);                 % and exit
            y = y + V(:,1:j)*z;
            % error(i+1) = abs(s(i+1)) / bnrm2;
            break;
        end
    end
    
    if ( error(itt+1) <= epsi ), break, end
    z = H(1:m,1:m)\s(1:m);
%        x = x + V*y;                                   % update approximation
    y = y + V(:,1:m)*z;                                 % update approximation
    % compute residual
    r = b - A*M*y;
    s(j+1) = norm(r);
    error(itt+1) = s(j+1) / normb;                        % check convergence
    if ( error(itt+1) <= epsi )
        break;
    end
    iter = iter+1;
end
iter = itt;
x = M*y;