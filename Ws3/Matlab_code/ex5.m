clc
close all
clear all
%% ---------------------- Matrix in ws1 ------------------------------- %%
% N = 100;
% beta = 0.9;
% alp = 2;
% A = createMatrix(N,beta,alp);
%% -------------------- Matrix laplace 1d ----------------------------- %%
% N = 100;
% A = lap1d(N);

%% -------------------- Matrix laplace 2d ----------------------------- %%
% n = 10; N = n*n;
% A = lap2d(n,n);

%% -------------------- Matrix market ----------------------------- %%
% ff = 'Matrices/hor__131.mtx';
% ff = 'Matrices/pde225.mtx';
% ff = 'Matrices/saylr4.mtx';
% ff = 'Matrices/sherman4.mtx';
ff = 'Matrices/tub100.mtx';
% [rows, cols, entries, rep, field, symm] = mminfo(ff) ;
[A,rows,cols,entries,rep,field,symm] = mmread(ff);
N = rows;



maxiter = 100;
tol = 10^(-14);
xex= 100*(rand(N,1)-0.5);
b = A*xex;

x0 = zeros(N,1);
%% %%%%========================== Exercise 2 ==========================%%%%

% Full GMRES method
m = N;
tic
[x1,res1,iter1] = GMRES(A,b,x0,m,maxiter,tol);
time1 = toc;
err1 = norm(x1 - xex) ;
% disp(['error GMRES = ', num2str(err1), ...
%     ', iter = ', num2str(iter1), ', time = ', num2str(time1)])

% Semi- GMRES method
m = 50;
sm = num2str(m) ;
tic
[x2,res2,iter2] = GMRES(A,b,x0,m,maxiter,tol);
time2 = toc;
err2 = norm(x2 - xex) ;
% disp(['error GMRES(' sm,  ') = ', num2str(err2), ...
%     ', iter = ', num2str(iter2), ', time = ', num2str(time2)])

% GMRES diag precondtioner
m = N;
M = diag(diag(A));
tic
[x3,res3,iter3] = Right_PRECGMRES( A,b,x0,M,m,maxiter,tol);
time3  = toc;
err3 = norm(x3 - xex) ;
% disp(['error GMRES diag preconditioner = ', num2str(err3),...
%     ', iter = ', num2str(iter3), ', time = ', num2str(time3)])

% GMRES LU(0) precondtioner
m = N;
[L,U] = ilu(A);
M = inv(L*U);
tic
[x4,res4,iter4] = Right_PRECGMRES( A,b,x0,M,m,maxiter,tol);
time4  = toc;
err4 = norm(x4 - xex) ;
% disp(['error GMRES LU(0) precondtioner = ', num2str(err4),...
%     ', iter = ', num2str(iter4), ', time = ', num2str(time4)])

% DEFLGMRESR
m = N; %l = N+1;
[L,U] = ilu(A);
M = inv(L*U);
k = 5; kmax = 15;
tic
% [x5,res5,iter5] = DEFLGMRESR(A,b,x0,M,m,l,maxiter,tol);
[x5,res5,iter5] = PrecDeflGMRES(A,b,x0,m,k,kmax,maxiter,tol);
time5 = toc;
err5 = norm(x5 - xex) ;
% disp(['error DEFLGMRESR = ', num2str(err5),...
%     ', iter = ', num2str(iter5), ', time = ', num2str(time5)])

s1 = 'Full GMRES';
s2 = ['GMRES(' sm,  ')'];
s3 = 'GMRES diag precondtioner';
s4 = 'GMRES LU(0) precondtioner';
s5 = ['DEFLGMRESR(' num2str(k) ', '  num2str(kmax) ')'];

semilogy(res1,'-p ', 'DisplayName',s1)
hold on
semilogy(res2(1:500),'--d ', 'DisplayName',s2)
semilogy(res3,'-^ ', 'DisplayName',s3)
semilogy(res4,'-o ', 'DisplayName',s4)
semilogy(res5,'-* ', 'DisplayName',s5)


legend(s1,s2,s3,s4,s5)

iteration = [iter1;iter2;iter3;iter4;iter5];
Err = [err1;err2;err3;err4;err5];
Time = [time1;time2;time3;time4;time5];
Name = {s1;s2;s3;s4;s5};
table(iteration,Err,Time,'RowNames',Name)
