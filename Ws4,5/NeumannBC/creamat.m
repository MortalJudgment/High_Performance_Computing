function [A] = creamat(M,dem)

% A Laplace Matrix + Dirichlet boundary conditions

m1 = M.Nx ;
m2 = M.Ny ;
mm = m1*m2 ;

e = ones(mm,1);
udx2 = 1./(M.hx*M.hx) ;
udy2 = 1./(M.hy*M.hy) ;
c2 = e*udy2 ;
c1 = e*udx2 ;
c0 = 2*(c1 + c2);
A = spdiags([-c2 -c1 c0 -c1 -c2], ...
    [-m1 -1 0 1 m1 ], mm, mm) ;

for i = 1:M.Ndir
    ii = M.ldir(i) ;
    A(ii,1:mm) = 0. ;
    A(ii,ii) = 1. ;
end


%% Add loop to take into account Neumann interface conditions

if dem == 1 % Backward
    for i = 1:M.Nintr
        ii = M.lintr(i);
        A(ii,1:mm) = 0. ;
        A(ii,ii) = 3./(2*M.hx) ;
        A(ii,ii-1) = -4./(2*M.hx) ;
        A(ii,ii-2) = 1./(2*M.hx) ;
    end
end
if dem == 2 % Forward
    for i = 1:M.Nintr
        ii = M.lintr(i);
        A(ii,1:mm) = 0. ;
        A(ii,ii) = -3./(2*M.hx) ;
        A(ii,ii+1) = 4./(2*M.hx) ;
        A(ii,ii+2) = -1./(2*M.hx) ;
    end
end

