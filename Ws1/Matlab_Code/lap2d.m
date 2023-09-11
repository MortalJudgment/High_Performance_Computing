function A=lap2d(nx,ny)
% A=lap2d(nx,ny) matrix of -delta in 2d on a grid 
% of  nx  internal points in x and ny  internal points in y
% numbered by row. Uses the function kron of matlab
Dxx=lap1d(nx);
Dyy=lap1d(ny);

A=kron(speye(size(Dyy)),Dxx) + kron(Dyy,speye(size(Dxx)));
