% function [x,res,iter] = CGprecd(C,A,b,x0,tol,maxiter)
% % Gradient conjugue preconditioning: 
% % [X,res,iter] = CGprecd(C,A,b,x0,tol,maxiter) computes an approximation
% % for the solution of the linear system Ax = b performing n steps 
% % of the conjugate gradient method starting with x0, 
% % and returns in res the norm of the residual
% normb = norm(b);
% if  ( normb == 0.0 )
%     normb = 1.0;
% end
% x = x0; 
% r = b - A*x;
% res(1) = norm(r)/normb;
% it = 1;
% z = C\r;
% p = z;
% while res(it)>tol  && it<maxiter
%   Ap = A*p;                   % Compute Ap^m
%   pAp = p'*Ap;                % inner product each iteration 
%   num = z'*r;
%   alpha = num/pAp;
%   x = x0 + alpha*p;
%   x0 = x;                      % update the initial value
%   r = r - alpha*Ap;
%   z = C\r;
%   beta = z'*r/num;
%   p = z + beta*p;
%   res(it+1)= norm(r)/normb;
%   it = it+1;
% end
% iter = it
function [x,res,iter]=CGprecd(C,A,b,x0,tol,maxiter)

%
% Gradient conjugue preconditionne: 
%                    [x,e]=CGprecd(A,b,x0,tol,maxiter) resout Ax=b en 
%                    utilisant l'algorithme du gradient conjugue 
%                    preconditionne ar une matrice C^{-1}.
%                    tol est la tolerance, x0 le vecteur initial, maxiter le
%                    nombre d'iteratios maximal, et e contient l'erreur 
%                    relative a chaque iteration.
%                    la matrice A doit etre symetrique definie positive.
%


r=b-A*x0;
res(1)=norm(r);
z=C\r;
p=z;
k=1;

while res(k)>tol  && k<maxiter

  Ap=A*p;                   % pour avoir seulement un produit 
  pAp=p'*Ap;                % matrice vecteur par iteration    
  num=z'*r;
  alpha=num/pAp;
  
  x=x0+alpha*p;
  x0=x;
  
  r=r-alpha*Ap;
  z=C\r;
  beta=z'*r/num;
  p=z+beta*p;
  res(k+1)= norm(r);
  k=k+1;

end

iter = k;