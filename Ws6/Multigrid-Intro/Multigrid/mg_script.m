% test script for multigrid solver
%
%  VERSION 1.0
%  DATE 25.3.2004
%  EMAIL bernd@flemisch.net

steps = 3;        % number of refinement steps
Hmax = 0.5;       % maximum initial element diameter
smooth = 'GS';    % type of smoother, 
                  % 'Jacobi' for damped Jacobi or 'GS' for Gauss-Seidel
omega = 0.7;      % relaxation parameter for the smoother
nu1 = 6;          % number of pre-smoothing steps
nu2 = 6;          % number of post-smoothing steps
mu = 1;           % number of recursive multigrid step calls, 
                  % 1 for V-cycle, 2 for W-cycle
tol = 1e-4;       % tolerance for stopping criterion
maxit = 100;      % maximum number of multigrid iterations
exact = 1;        % flag whether exact solution is available
direct = 1;       % flag whether a direct solve is performed
prec = 0;         % flag whether multigrid is used as preconditioner
geom = 'squareg'; % Matlab geometry description
b = 'squareb1';   % Matlab boundary condition description
a = 1;            % problem parameters for the PDE toolbox to solve
c = 0;            % div(a*grad(u)) + c*u = f
f = '5*pi*pi*sin(2*pi*x).*sin(pi*y)'; % right hand side
u = 'sin(2*pi*x).*sin(pi*y)'; % exact solution
plot_fl = 1;      % flag whether to plot the FE solution

% [p1, e1, t1] = initmesh(geom, 'Hmax', Hmax);
% pdemesh(p1,e1,t1) 

[p, e, t, u_lmax , matrices] = mg_main(steps, Hmax, smooth, omega, nu1, nu2, mu, tol, maxit, ...
			    exact, direct, prec, geom, b, a, c, f, u);
% figure(2)
% pdemesh(p,e,t) 

if (plot_fl) 
  pdesurf(p, t, u_lmax);
  colormap('jet');
end