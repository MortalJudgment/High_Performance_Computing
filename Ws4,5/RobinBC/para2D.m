%Seq2D
clear all
close all
clc
global drecs2

%Split Topology 
T.Ndomx = 3;
T.Ndomy = 1;
T.Ndom = T.Ndomx * T.Ndomy ;
Ndom = T.Ndom ;

%----------------------------
%Initialisation des structures
%----------------------------

for i=1:Ndom
    M(i) = struct('Nx',0,'Ny',0,'a1',0,'a2',0,'b1',0,'b2',0, ...
        'hx',0,'x',0,'hy',0,'y',0, ...  % Gometry
        'Ndom',Ndom,'idom',i,  ...      % Number of dom. and  Domain number.
        'Ndir',0,'ldir',0,  ...         % Number Dir. points and list
        'Ninte',0,'linte',0, ...        % Number of sending Interface points  and list
        'Nintr',0,'lintr',0, ...        % Number de receiving Interface points and list
        'nxe',0,'nye',0,...             % Inward Normal for sending list
        'nxr',0,'nyr',0,...             % Outward Normal for receiving  list
        'Nin',0,'in',0, ...             % Number and list of inner points
        'Nvois',0,'novois',0,...        % Number of neighbours, and list 
        'lvoise',0,'lvoisr',0,...       % Concatenated  List of send and receive points
        'list',0);                      % Position in the list 
                                        % Exemple 3 sub domain split with  3 interface points
                                        % Domain 2 has 2 neighbouring domains (3,1) 
                                        % data to send to domain 1 are in second position
                                        % data to send to domain 3 are in first position


    P(i) = struct('usol',0,'uana',0,'udir',0,'sm0',0,'sm',0, ...
        'A',0) ;
end
%----------------------------
%Parameters
%----------------------------
%Parameters  global mesh
for n = 2:1:9
drecs2 = n ;   % .5*Number of overlap cells

Mg.a1 = 0;    %  global box
Mg.a2 = 2;
Mg.b1 = 0;
Mg.b2 = 2;

Mg.Nx = 41 ;  % Total number of points
Mg.Ny = 21 ;

% Subdomain Indices
for i = 1:T.Ndom
    T.l(i,1) = mod(i-1,T.Ndomx) + 1 ;
    T.l(i,2) = (i- T.l(i,1))/T.Ndomx + 1 ;
    % Only a split in the x direction is performed in this version
%     T.l(i,1) = i ;
%     T.l(i,2) = 1 ;
end

%----------------------------
%Geometry
%----------------------------
for i=1:Ndom
    [M(i), T] = mcreamesh(Mg, M(i), T, i) ;
end
% visu(1,T,M,P) ;

%----------------------------
% matrix , rhs
%----------------------------
%% Robin Boundary condition
% x_u^n_1 + rho*u^n_1 = x_u^n_2 + rho*u^n_2

rho = 1;

for i=1:Ndom
    % rhs  +  Dirichlet boundary conditions
    [P(i).sm0, P(i).udir, P(i).uana, P(i).usol] = init(M(i), Mg, 1) ;
    %Matrix
%     [P(i).A] = creamat(M(i),i,rho) ;
    [P(i).A] = creamatrix(M(i),rho);
end

%----------------------------
% Solve
% %----------------------------
maxit= 500 ;
tol = 1.e-6 ;
% P = schwarz(P,M,T,Mg,tol,maxit);
[P,res] = schwarz_Jacobi(P,M,T,rho,Mg,tol,maxit);
%----------------------------
% Error
%----------------------------
umax = -1. ;
u2 = 0. ;
for i=1:Ndom
    P(i).uu = P(i).uana - P(i).usol;
    umax = max( umax, max(abs(P(i).uu))) ;
    u2 = u2 + dot(P(i).uu,P(i).uu);
end
er2 = sqrt(u2/(Mg.Nx*Mg.Ny)) ;
sprintf('er2= %e  ,erinf= %e\n', er2, umax)
num = num2str(n);
semilogy(res,'-d','DisplayName',num)
hold on
end

lgd = legend('show');
title(lgd,'Number of overlap cells') 
visu(2,T,M,P);
visu(3,T,M,P);



