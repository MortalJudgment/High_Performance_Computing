function [sm0, udir, uana, usol] = init(M,Mg,icase)

m1 = M.Nx ;
m2 = M.Ny ;
mm = M.Nx*M.Ny ;
x = M.x ;
y = M.y ;
[xx,yy] = meshgrid(M.x,M.y) ;
xx = xx' ;
yy = yy' ;

switch icase
    case 1  % u = sin(om*x + cx)*sin(om*y + cy)+ x*x + y*y
        omx = 2*pi/(Mg.a2 - Mg.a1) ;
        cx = -omx * Mg.a1 ;
        ux = sin(omx * x' + cx) ;

        omy = 2*pi/(Mg.b2 - Mg.b1) ;
        cy = -omy * Mg.b1 ;
        uy = sin(omy * y' + cy) ;

        uu = ux *uy' ;
        uana = reshape(uu + xx.*xx + yy.* yy,mm,1) ;
        sm0 = reshape( (omx*omx+ omy*omy)*uu - 4 ,mm,1) ;
        
    case 2  % u = x + y ;
        uana = reshape(xx + yy,mm,1)  ;
        sm0 = zeros(mm,1) ;
    case 22  % u = x ;
        uana = reshape(xx,mm,1)  ;
        sm0 = zeros(mm,1) ;
                
    case 3  % u = x*x + y*y ;
        uana = reshape(xx.*xx + yy.* yy,mm,1)  ;
        sm0 = -4.*ones(mm,1) ;

    otherwise
        disp(['unknown case' ])
        sm0 = 0 ;
        uana= 0 ;
        udir = 0;
end
udir = zeros(mm,1) ;

if(M.Ndir ~= 0)
    udir(M.ldir) = uana(M.ldir) ;
    sm0(M.ldir) = uana(M.ldir) ;    % apply Dirichlet boundary condition to sm0 
end
usol = zeros(mm,1) ;
