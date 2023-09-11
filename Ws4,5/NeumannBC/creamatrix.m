function [A] = creamatrix(M)

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

% Boundary condition
for i = 1:M.Ndir
    ii = M.ldir(i) ;
    A(ii,1:mm) = 0. ;
    A(ii,ii) = 1. ;
end


%% Add loop to take into account Neumann interface conditions
% list of interface boundary points to be updated
for j = 1:M.Nvois
    indSi = (M.lvoisr(j)+1 : M.lvoisr(j+1));
    indAi = M.lintr(indSi) ;
    n = length(indAi);
    
    for k = 1:n
        ii = indAi(k);
        % Vector face take value -1 or 1
        point = M.nxr((j-1)*n + k);
        A(ii,1:mm) = 0. ;
        A(ii,ii)            = point*3./(2*M.hx) ;
        A(ii,ii - point*1)  = point*(-4)./(2*M.hx) ;
        A(ii,ii - point*2)  = point*1./(2*M.hx) ;
%         A(ii,:)
    end

end



