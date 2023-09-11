clc
close all
clear all
%%
n = 10;
N = n*n;
% A = lap1d(N);
A = lap2d(n,n);
A = A'*A;
x_ex = rand(N,1);
b = A*x_ex;


x0 = zeros(N,1);
tol = 10^(-8);
% tol = 10^(-14);
itmax = 400;
%% %%%%========================== Exercise 1 ==========================%%%%%
% Jacobi method
tic
[x1,res1,it1,time1] = Jacobi(A,b,x0,tol,itmax);
time = toc;
err = norm(x1 - x_ex ) ;
disp(['error Jacobi = ', num2str(err), ', iter = ', num2str(it1), ', time = ', num2str(time)])

% GaussSeidel
tic 
[x2,res2,it2,time2] = GaussSeidel(A,b,x0,tol,itmax);
time = toc;
err = norm(x2 - x_ex ) ;
disp(['error GS = ', num2str(err), ', iter = ', num2str(it2), ', time = ', num2str(time)])

% Backward GaussSeidel
tic
[x3,res3,it3,time3] = Backward_GaussSeidel(A,b,x0,tol,itmax);
time = toc;
err = norm(x3 - x_ex ) ;
disp(['error BGS = ', num2str(err), ', iter = ', num2str(it3), ', time = ', num2str(time)])

% Simple projection method : Steepest Descent
tic
[x4,res4,it4,time4] = Steepest_Descent(A,b,x0,tol,itmax);
time = toc;
err = norm(x4 - x_ex ) ;
disp(['error SD = ', num2str(err), ', iter = ', num2str(it4), ', time = ', num2str(time)])

% Successive Over Relaxation (SOR)
i = 1.05;
tic
[x5,res5,it5,time5] = SOR(A,b,i,x0,tol,itmax);
time = toc;
err = norm(x5 - x_ex ) ;
disp(['error SOR = ', num2str(err), ', iter = ', num2str(it5), ', time = ', num2str(time)])

% Conjugate Gradient method (CG)
tic
[x6,res6,it6,time6] = CG(A,b,x0,tol,itmax);
time = toc;
err = norm(x6 - x_ex ) ;
disp(['error CG = ', num2str(err), ', iter = ', num2str(it6), ', time = ', num2str(time)])

% s1 = 'Jacobi' ;
% semilogy(res1,'-*', 'DisplayName',s1)
% hold on
s2 = 'Gauss-Seidel' ;
semilogy(res2,'-*', 'DisplayName',s2)
hold on
s3 = 'Backward-Gauss-Seidel' ;
semilogy(res3,'-*', 'DisplayName',s3)
s4 = ' Steepest Descent' ;
semilogy(res4,'-*', 'DisplayName',s4)
s5 = ' SOR' ;
semilogy(res5,'-*', 'DisplayName',s5)
s6 = ' CG' ;
semilogy(res6,'-*', 'DisplayName',s6)


legend(s2,s3,s4,s5,s6)
title('Residual of different methods')
