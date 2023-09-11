clc
close all
clear all
%% -------------------------------------------------------------------- %%
% Matrix in ws1
% N = 100;
% beta = 0.9;
% alp = 2;
% A = createMatrix(N,beta,alp);

% Matrix of Laplace 1D equation
% N = 100;
% A = lap1d(N);

% Matrix of Laplace 2D equation
% n = 10; N = n*n;
% A = lap2d(n,n);

% ff = 'Matrices/hor__131.mtx';
% ff = 'Matrices/pde225.mtx';
% ff = 'Matrices/saylr4.mtx';
% ff = 'Matrices/sherman4.mtx';
ff = 'Matrices/tub100.mtx';
% [rows, cols, entries, rep, field, symm] = mminfo(ff) ;
[A,rows,cols,entries,rep,field,symm] = mmread(ff);
N = rows;

A = A'*A ;

xex= rand(N,1);
b = A*xex;

x0 = zeros(N,1);
maxiter = 1000;
tol = 10^(-12);
% %%%%%========================== Exercise 1 ==========================%%%%%
% % Exercice 1 - Jacobi
% [u_Jac,error0,iter]=Jacobi(A,b,x0,maxiter,tol);
% err = norm(u_Jac - U_ex ) ;
% disp(['error Jacobi = ', num2str(err), ', iter = ', num2str(iter)])
% 
% % Exercice 1 - Gauss-Seidel lower
tic
[u_gsl,error1,iter] = GSL(A,b,x0,maxiter,tol);
time = toc;
err = norm(u_gsl - xex ) ;
disp(['error GSL = ', num2str(err), ', iter = ', num2str(iter), ', time = ', num2str(time)])

% Exercice 1 - Gauss-Seidel upper
[u_gsu,error2,iter] = GSU(A,b,x0,maxiter,tol);
err = norm(u_gsu - xex ) ;
disp(['error GSU = ', num2str(err), ', iter = ', num2str(iter), ', time = ', num2str(time)])

% Exercice 1 - Steepest Descent
tic
[u_std,error3,iter] = STD(A,b,x0,maxiter,tol);
time = toc;
err = norm(u_std - xex ) ;
disp(['error STD = ', num2str(err), ', iter = ', num2str(iter), ', time = ', num2str(time)])

% Exercice 1 - SOR
omg = 1.2 ;
tic
[u_sor,error4,iter] = SOR(A,b,x0,maxiter,tol, omg);
time = toc;
err = norm(u_sor - xex ) ;
disp(['error SOR = ', num2str(err), ', iter = ', num2str(iter), ', time = ', num2str(time)])

% Exercice 1 - Conjugate Gradient
tic
[u_cg,error5,iter] = CG(A,b,x0,maxiter,tol);
time = toc;
err = norm(u_cg - xex ) ;
disp(['error CG = ', num2str(err), ', iter = ', num2str(iter), ', time = ', num2str(time)])

% Exercice 1 - GMRES
tic
[u_gmres,error6,iter] = GMRES(A,b,x0,N,maxiter,tol) ;
time = toc;
err = norm(u_gmres - xex ) ;
disp(['error GMRES = ', num2str(err), ', iter = ', num2str(iter), ', time = ', num2str(time)])


% Visualize convergence versus iterations
% s0 = 'Jacobi' ;
% semilogy(error0,'-x','DisplayName',s0)
s1 = 'Gauss-Seidel Lower' ;
semilogy(error1,'-x', 'DisplayName',s1)
hold on
s2 = 'Gauss-Seidel Upper' ;
semilogy(error2,'-x', 'DisplayName',s2)
s3 = ' Steepest Descent' ;
semilogy(error3,'-x', 'DisplayName',s3)
s4 = ' SOR' ;
semilogy(error4,'-x', 'DisplayName',s4)
s5 = ' CG' ;
semilogy(error5,'-x', 'DisplayName',s5)
s6 = ' GMRES' ;
semilogy(error6,'-x', 'DisplayName',s6)
legend(s1,s2,s3,s4,s5,s6)



