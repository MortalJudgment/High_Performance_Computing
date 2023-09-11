clc
close all
clear all
%%
alp = 2;
beta = 0.9;
n = 100;
A = create_matrix(n,beta,alp)  ;

x_ex = rand(100,1);
b = A*x_ex;

x0 = zeros(100,1);
tol = 10^(-12);
itmax = 100;
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
i = 1.2;
tic
[x5,res5,it5,time5] = SOR(A,b,i,x0,tol,itmax);
time = toc;
err = norm(x5 - x_ex ) ;
disp(['error SOR = ', num2str(err), ', iter = ', num2str(it5), ', time = ', num2str(time)])

s1 = 'Jacobi' ;
semilogy(res1,'-*', 'DisplayName',s1)
hold on
s2 = 'Gauss-Seidel' ;
semilogy(res2,'-*', 'DisplayName',s2)
s3 = 'Backward-Gauss-Seidel' ;
semilogy(res3,'-*', 'DisplayName',s3)
s4 = ' Steepest Descent' ;
semilogy(res4,'-*', 'DisplayName',s4)
s5 = ' SOR' ;
semilogy(res5,'-*', 'DisplayName',s5)

legend(s1,s2,s3,s4,s5)
title('Residual of different methods')





%% Compare different value of relax in SOR method
% % Successive Over Relaxation (SOR)
close all
er = [];
iter = [];
timeC = [];
c = 0:0.05:2;
for ii = c
    tic
    [x,res,it] = SOR(A,b,ii,x0,tol,itmax);
    time = toc;
    err = norm(x - x_ex ) ;
    er = [er;err];
    iter = [iter;it];
    timeC = [timeC;time];
    semilogy(res,'-*')
    hold on
end
Relax = c';
table(Relax,er,iter,timeC)
