%************************************************
function[x,error,iter]=PrecDeflGMRES(A,b,x,m,k,kmax,itmax,epsi)
%************************************************

% gmres.m solves the linear system Ax=b
% using the Generalized Minimal residual ( GMRESm ) method with restarts .
%
% input   M        REAL nonsymmetric positive definite preconditionning matrix
%         A        REAL nonsymmetric positive definite matrix
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

normb = norm( b );
if  ( normb == 0.0 )
    normb = 1.0;
end
n= size(A,1);
e=ones(n,1) ;
M = spdiags(e, 0, n,n) ;
invT = [] ;

% residual
r=b-A*x;
error(1) = norm( r )/normb;

if ( error(1) < epsi )
    return;
end

V(1:n,1:m+1) = zeros(n,m+1);
H(1:m+1,1:m) = zeros(m+1,m);   % Hessenberg matrix
Hess(1:m,1:m) = zeros(m,m);   % Hessenberg matrix
cs(1:m) = zeros(m,1);
sn(1:m) = zeros(m,1);
e1    = zeros(n,1); % basic vector
e1(1) = 1.0;
iter=1;  % step of iterator
itt = 0 ;
U=[];
T=[];
tt=0. ;
while iter <= itmax                              % begin iteration
    r=b-A*x;
    V(:,1)=r/norm(r);
    s = norm(r)*e1;
    for j = 1:m,
        itt = itt+ 1 ;% construct orthonormal
        % basis using Gram-Schmidt
        %         w=A*(M*V(:,j));
        w=A*Mprod(iter,V(:,j),invT,U )  ;
        for i = 1:j
            H(i,j)= w'*V(:,i);
            w = w - H(i,j)*V(:,i);
        end
        H(j+1,j) = norm(w);
        V(:,j+1) = w/H(j+1,j);
        Hess(1:j+1,j) = H(1:j+1,j) ;
        % We tranform the Hessenberg matrix H into a triangular matrix by applying Givens rotation
        for i = 1:j-1,
            temp     =  cs(i)*H(i,j) + sn(i)*H(i+1,j);
            H(i+1,j) = -sn(i)*H(i,j) + cs(i)*H(i+1,j);
            H(i,j)   = temp;
        end
        [cs(j),sn(j)] = rotmat( H(j,j), H(j+1,j) ); % form j-th rotation matrix
        temp   = cs(j)*s(j);
        s(j+1) = -sn(j)*s(j);
        s(j) = temp;
        H(j,j) = cs(j)*H(j,j) + sn(j)*H(j+1,j);
        H(j+1,j) = 0.0;  %eliminate H(j+1,j)
        error(itt+1) = abs(s(j+1)) / normb;
        if ( error(itt+1) <= epsi ),                        % update approximation
            y = H(1:j,1:j) \ s(1:j);                 % and exit
            %             x = x + M*(V(:,1:j)*y);
            yy =V(:,1:j)*y ;
            x = x + Mprod(iter, yy,invT,U )  ;
            % error(i+1) = abs(s(i+1)) / bnrm2;
            break;
        end
    end
    
    if ( error(itt+1) <= epsi )
        break;
    end
    y = H(1:m,1:m)\s(1:m);
    %    x = x + V*y;                                   % update approximation
    %     x = x + M*(V(:,1:j)*y);                                 % update approximation
    yy =V(:,1:j)*y ;
    x = x + Mprod(iter,yy,invT,U )  ;
    
    % compute residual
    r=b-A*x;
    s(itt+1) = norm(r);
    error(itt+1) = s(itt+1) / normb;                        % check convergence
    if ( error(itt+1) <= epsi )
        break;
    end
    iter=iter+1;
    
    %     tic ;
    % Calcul de la nouvelle matrice de deflation
    [E,R] = schur(Hess(1:m,1:m));
    I=sort(diag(R));
    X=[];
    
    for i=1:k
        for j=1:m
            if (I(i)==R(j,j))
                X=[X E(:,j)];
                break;
            end
        end
    end
    W=V(:,1:m)*X;
    
    rmax0 = size(U,2) ;
    if(rmax0 <= kmax)
        if(isempty(U)==1)   % % % %Si U vide
            U1 = W ;
            AU1 = A*U1 ;
            T=U1'*AU1;
            rmax=size(U1,2);
            U=U1 ;
        else
            U= U1 ;
            for i=1:k
                for j = 1:size(U,2)
                    dd= W(:,i)'*U(:,j);
                    W(:,i) = W(:,i) - dd*U(:,j);
                end
                W(:,i) = W(:,i)/norm(W(:,i)) ;
                U = [U W(:,i)] ;
            end
            rmax=size(U,2);
            U2 = U(:,rmax0+1:rmax) ;
            AU2 = A*U2 ;
            
            T = [T U1'*AU2] ;
            U1 = U ;
            AU1 = [AU1 AU2] ;
            T(rmax0+1:rmax,1:rmax) = U2'*AU1 ;
        end
    end
    
    invT = abs(I(m))*inv(T) ;
    %         M=eye(n,n)+U*(abs(I(m))*inv(T)-eye(rmax,rmax))*U';
    %      tt = tt + toc ;
end

% x = M*x ;
iter = itt ;
%     tt