function [M,T] = mcreamesh(Mg, M, T, i)

global drecs2            % Number of overlap cells
drec = 2* drecs2 ;
drec1 = drec + 1 ;

% Splitting Topology
%---------------------------------------------
Ndom = T.Ndom;
Ndomx = T.Ndomx ;
Ndomy = T.Ndomy ;
ix = T.l(i,1)  ;
iy = T.l(i,2) ;


% Number of points per domain and step size hx and hy
%---------------------------------------------------
npmx = floor((Mg.Nx - 1)/Ndomx) ;  % Number of points per domain (as balanced as possible)
nrx = Mg.Nx -1 - npmx*Ndomx;
M.Nx = npmx + 1 + (ix <=nrx) ;

npmy = floor((Mg.Ny - 1)/Ndomy) ;   % Number of points per domain (as balanced as possible)
nry = Mg.Ny -1 - npmy*Ndomy;
M.Ny = npmy + 1 + (iy <=nry) ;

hx = (Mg.a2 - Mg.a1)/(Mg.Nx-1) ;   % Common stepsize
% otherwise problem with normal vector
hy = (Mg.b2 - Mg.b1)/(Mg.Ny-1) ;   % Common stepsize
% otherwise problem with normal vector

%Split  [a1,a2] x [b1,b2] with no overlap
%---------------------------------------------------
M.a1 = Mg.a1 ;
i0=0 ;
for j=1:ix-1
    M.a1 = M.a1 + hx*( npmx  + (j <=nrx) ) ;
    i0=i0 + npmx  + (j <=nrx) ;
end
M.a2 = M.a1 + hx*(M.Nx-1) ;

M.b1 = Mg.b1 ;
j0 = 0 ;
for j=1:iy-1
    M.b1 = M.b1 + hy*( npmy  + (j <=nry) ) ;
    j0=j0 + npmy  + (j <=nry) ;
end
M.b2 = M.b1 + hy*(M.Ny-1) ;

T.i0 = i0 ;
T.j0 = j0 ;

% Overlap
%------------------
if  (Ndom > 1)
    if(ix > 1)
        M.a1 = M.a1-drecs2*hx ;
        M.Nx = M.Nx +drecs2 ;
    end
    if(ix < Ndomx)
        M.a2 = M.a2+drecs2*hx;
        M.Nx = M.Nx +drecs2 ;
    end
end


%Mesh
M.hx = hx ;
M.x(1:M.Nx) = M.a1 + (0:1:M.Nx-1)*M.hx ;

M.hy = hy ;
M.y(1:M.Ny) = M.b1 + (0:1:M.Ny-1)*M.hy ;
                             %   l13   l14      l15      l16
M.Ndom = Ndom ;              %     ------------------------
M.idom = i ;                 %    |     |        |        |
mm = M.Nx*M.Ny ;             %    |     |        |        |
mx = M.Nx;                   %    |     |        |        |
my = M.Ny ;                  %    |     |        |        |
icx = 1;                     %    |     |        |        |
icy = mx ;                   %    |     |        |        |
l1 = 1 ;                     %    |     |        |        |
l2 = drec1 ;                 %    |     |        |        |
l3 = mx-drec ;               %    |     |        |        |
l4 = mx ;                    %    |     |        |        |
                             %     ------------------------
                             %   l1    l2       l3       l4
l13 = mm - mx+1 ;
l14 = l13 + drec ;
l15 = mm - drec ;
l16 = mm ;

% % Dirichlet points and interface points
M.Ndir = 0 ;
M.Nvois = 0 ;
M.Nintr = 0 ;
M.lvoisr(1) = 0 ;
M.Ninte = 0 ;
M.lvoise(1) = 0 ;
if(M.Ndom == 1)
    % Dirichlet points
    M.Ndir = 2*mx + 2*(my-2);
    M.ldir(1:mx) = (l1 : icx : l4) ;
    M.ldir(mx + 1: mx + my - 2) = (l4+icy : icy : l16-icy) ;
    M.ldir( mx + my - 1: 2*mx + my-2 ) = ( l13 : icx : l16) ;
    M.ldir( 2*mx + my - 1 : 2*mx + 2*(my-2)) = (l1+icy : icy : l13-icy) ;
else
    % % Dirichlet points and interface points
    if(iy == 1) % bottom points
        M.ldir(M.Ndir+1:M.Ndir+mx) = (l1: icx :l4) ;
        M.Ndir = M.Ndir + mx ;
    end
    if(ix == Ndomx) % bottom points
        M.ldir(M.Ndir+1 : M.Ndir + my-2) = (l4 + icy : icy : l16-icy) ;
        M.Ndir = M.Ndir + my-2 ;
    end
    if(iy == Ndomy)
        M.ldir(M.Ndir+1: M.Ndir + mx) = ( l13: icx :l16) ;
        M.Ndir = M.Ndir + mx ;
    end
    if(ix == 1)
        M.ldir(M.Ndir+1 : M.Ndir + my-2) = (l1+icy : icy : l13-icy) ;
        M.Ndir = M.Ndir + my-2 ;
    end
    % Overlap Points send and receive
    %---------   Face South -----------------------------------
    % DIRICHLET
    
    %---------  Fin  Face South -----------------------------------
    
    %---------  Face East -----------------------------------
    if(ix < Ndomx)
        M.Nvois = M.Nvois + 1 ;
        nv = (iy-1)*Ndomx + ix +1 ;
        M.novois(M.Nvois) = nv ;
        M.list(nv) = M.Nvois ;
        
        % reception l4,l16
        % send     l3,l15
        
        M = FaceInt(M,l4,l16,l3,l15,icy, -icx,1,0) ;
        
    end
    %---------  End  Face Est -----------------------------------
    
    %---------  Face Nord -----------------------------------
    % DIRICHLET
    %---------  End  Face Nord -----------------------------------
    
    %---------  Face West-----------------------------------
    if(ix > 1)
        M.Nvois = M.Nvois +1 ;
        nv = (iy-1)*Ndomx + ix -1 ;
        M.novois(M.Nvois) = nv ;
        M.list(  nv  ) = M.Nvois ;
        
        % reception l1,l13
        % Send     l2,l14
        
        M = FaceInt(M,l1,l13,l2,l14,icy,icx,-1,0) ;
        
    end
    
    %---------  End  Face Ouest -----------------------------------
    
end

% List of inner points
M.Nin = (mx-2)*(my-2);
Lin = ones(1,mx*my) ;
if(M.Nintr ~= 0)  Lin(M.lintr) = 0 ; end
if(M.Ndir ~= 0)  Lin(M.ldir) = 0 ; end
num = (1:1:mx*my)  ;
M.in = num(Lin==1) ;
%
