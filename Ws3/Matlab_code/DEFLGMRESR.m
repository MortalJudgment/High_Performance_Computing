%************************************************
function[x,error,iter, itv] = DEFLGMRESR(A,b,x,M,m,l,itmax,epsi)  
%************************************************
% DEFLGMRESR.m solves the linear system Ax=b
% using the Generalized Minimal residual ( GMRESm ) method with restarts .
%  With Right preconditioning and Schur decomposition
%
% input   A        REAL nonsymmetric positive definite matrix
%         x        REAL initial guess vector
%         bb        REAL right hand side vector
%         m        INTEGER number of iterations between restarts
%         l
%         itmax   INTEGER maximum number of iterations
%         epsi      REAL error tolerance
%         M         Right Preconditioner
%
% output  x        REAL solution vector
%         error    REAL error norm
%         iter     INTEGER number of iterations performed
%         itv Total number of iterations


% initialization
normb = norm( b );
if  ( normb == 0.0 ),
    normb = 1.0;
end
% invM = inv(M);
% residual
r = b - A*(M\x);
error(1) = norm( r )/normb;

if ( error(1) < epsi ) 
    return; 
end

n= size(A,1);                                  
V(1:n,1:m+1) = zeros(n,m+1);   
H(1:m+1,1:m) = zeros(m+1,m);   % Hessenberg matrix
U = [];
cs(1:m) = zeros(m,1);
sn(1:m) = zeros(m,1);
e1    = zeros(n,1); % basic vector
e1(1) = 1.0;
iter=1;  % step of iterator
itv= 0 ; % Total iterations
while iter <= itmax                             % begin iteration
      r = b - A*(M\x);
      V(:,1)=r/norm(r);
      s = norm(r)*e1;
   for j = 1:m                                  % construct orthonormal
         itv = itv +1 ;                         % basis using Gram-Schmidt
       w = A*(M\V(:,j));
	   for i = 1:j
	       H(i,j)= w'*V(:,i);
	       w = w - H(i,j)*V(:,i);
       end
	   H(j+1,j) = norm(w);
	   V(:,j+1) = w/H(j+1,j);
   % We tranform the Hessenberg matrix H into a triangular matrix by applying Givens rotation 
	   for i = 1:j-1                              
           temp     =  cs(i)*H(i,j) + sn(i)*H(i+1,j);
           H(i+1,j) = -sn(i)*H(i,j) + cs(i)*H(i+1,j);
           H(i,j)   = temp;
       end
	   [cs(j),sn(j)] = rotmat( H(j,j), H(j+1,j) ); % form i-th rotation matrix
       temp   = cs(j)*s(j);                        
       s(j+1) = -sn(j)*s(j);
	   s(j) = temp;
       H(j,j) = cs(j)*H(j,j) + sn(j)*H(j+1,j);
       H(j+1,j) = 0.0;  %eliminate H(j+1,j)
 	   error(j+1) = abs(s(j+1)) / normb;
 	   if ( error(j+1) <= epsi )                   % update approximation
	      y = H(1:j,1:j) \ s(1:j);                 % and exit
          x = x + V(:,1:j)*y;
        % error(i+1) = abs(s(i+1)) / bnrm2;
	      break;
	   end
   end
 
   if ( error(j+1) <= epsi ), break, end
   y = H(1:m,1:m)\s(1:m);
%    x = x + V*y;                                   % update approximation
   x = x + V(:,1:j)*y;                              % update approximation
   
   % Compute Shur vectors of H noted Sl, l = 2
   [eigvals] = eig(H(1:m,1:m));
   eigval_new = [];
   for i = 1:length(eigvals)
       if (eigvals(i) ~= min(eigvals))
           eigval_new = [eigval_new;eigvals(i)];
       end
   end
   min_eigval = [min(eigvals); min(eigval_new)];
   [U_Hess,T] = schur(H(1:m,1:m));
   Sl = zeros(m,l);
   for i = 1:l
       for jj = 1:size(T,2)
           if T(jj,jj) == min_eigval(i)
               Sl(:,i) = T(:,jj);
           end
       end    
   end
   k = size(U,2);
   W = V(:,1:m)*Sl;
   U = [U W];
   % Orthogonalize VmSl against U
   U(:,1) = U(:,1) / sqrt(U(:,1)'*U(:,1));
   for i = 2:k
       for ii = 1:k-1
           U(:,i) = U(:,i) - ( U(:,ii)'*U(:,i) )/( U(:,ii)'*U(:,ii) )*U(:,ii);
       end
       U(:,i) = U(:,i) / sqrt(U(:,i)'*U(:,i));
   end
  T = U'*A*U;
  invM = eye(n,n) + U*(abs(max(eig(full(A))))*inv(T) - eye(size(T,1)))*U';
   
   % compute residual
   r = b - A*(M\x);
   s(j+1) = norm(r);
   error(j+1) = s(j+1) / normb;                        % check convergence
   if ( error(j+1) <= epsi ) 
       break; 
   end
   iter=iter+1;
end
x = M\x;    % true solution after preconditioning
