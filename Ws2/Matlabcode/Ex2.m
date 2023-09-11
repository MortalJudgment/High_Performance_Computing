% %%%%%========================== Exercise 2 ==========================%%%%%
clc
close all
clear all

%% Matrix in ws1
% N = 100;
% beta = 0.9;
% alp = 2;
% A = createMatrix(N,beta,alp);
%% Matrix in market
% ff = 'Matrices/hor__131.mtx';
% ff = 'Matrices/pde225.mtx';
% ff = 'Matrices/saylr4.mtx';
% ff = 'Matrices/sherman4.mtx';
ff = 'Matrices/tub100.mtx';
% [rows, cols, entries, rep, field, symm] = mminfo(ff) ;
[A,rows,cols,entries,rep,field,symm] = mmread(ff);
N = rows;

xex= rand(N,1);
b = A*xex;

x0 = zeros(N,1);
maxiter = 1000;
tol = 10^(-15);

% GMRES method
tic
[x1,res1,iter1] = GMRES(A,b,x0,N,maxiter,tol);
time = toc;
err = norm(x1 - xex) ;
disp(['error GMRES = ', num2str(err), ', iter = ', num2str(iter1), ', time = ', num2str(time)])

% GMRES right precondtioning method
[L,U] = ilu(A);
M = inv(L*U);
% M = diag(diag(A));
tic
[x2,res2,iter2] = Right_PRECGMRES( A,b,x0,M,N,maxiter,tol);
time = toc;
err = norm(x2 - xex) ;
disp(['error GMRES right preconditioner = ', num2str(err),...
    ', iter = ', num2str(iter2), ', time = ', num2str(time)])

s1 = 'GMRES';
semilogy(res1,'--pb ', 'DisplayName',s1)
hold on
s2 = 'GMRES Right Preconditioner';
semilogy(res2,'--dr ', 'DisplayName',s2)
legend(s1,s2)

