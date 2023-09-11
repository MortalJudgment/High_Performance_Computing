% TP3
%% ------------------------------ %%
clc
close all
clear all
%% ========================== Matrix A in ws1 ==========================%%
% N = 100;
% beta = 0.9;
% alp = 2;
% A = createMatrix(N,beta,alp);
% A = A*A';
%% ====================== sparse symetric matrix ======================%%
% n = 50; N = n;
% A = lap1d(n);


%% ====================== sparse symetric matrix ======================%%
% n = 10; N = n*n;
% A = lap2d(n,n);

%% ====================== matrix in market ======================%%
% ff = 'Matrices/hor__131.mtx'
% ff = 'Matrices/pde225.mtx';
ff = 'Matrices/saylr4.mtx'
% ff = 'Matrices/sherman4.mtx'
% ff = 'Matrices/tub100.mtx'
% [rows, cols, entries, rep, field, symm] = mminfo(ff) ;
[A,rows,cols,entries,rep,field,symm] = mmread(ff);
N = rows;
A = A'*A ;
%% ========================== input argument =========================%%
xex = rand(N,1);
b = A*xex;
x0 = zeros(N,1);
tol = 10^(-12);
maxiter = N;
%% ---------------------------------------------------------------------%%
% Conjugate gradient
In = eye(N,N);
[xCG,res0,iter0] = CGprecd(In,A,b,x0,tol,maxiter);
err = norm(xCG - xex) ;
disp(['error CG = ', num2str(err), ', iter = ', num2str(iter0)])

% Preconditioner CG with Gauss Seidel
om = 1; % Gauss Seidel
E = -tril(A);
vecdiag = diag(A); rootvec = 1./sqrt(vecdiag);
D = diag(vecdiag);
RD = diag(rootvec);

S = (D - om*E)*RD/sqrt(om*(2 - om));
CGS = S*S';
[x1,res1,iter1] = CGprecd(CGS,A,b,x0,tol,maxiter);
err = norm(x1 - xex) ;
disp(['error preconditioner CG with Gauss Seidel = ', num2str(err), ', iter = ', num2str(iter1)])

% % Preconditioner CG with SSOR
om = 1.2;
E = -tril(A);
vecdiag = diag(A); rootvec = 1./sqrt(vecdiag);
D = diag(vecdiag);
RD = diag(rootvec);

S = (D - om*E)*RD/sqrt(om*(2 - om));
CSSOR = S*S';
[x2,res2,iter2] = CGprecd(CSSOR,A,b,x0,tol,maxiter);
err = norm(x2 - xex) ;
disp(['error preconditioner CG with SSOR = ', num2str(err), ', iter = ', num2str(iter2)])

% % Preconditioner CG with Cholewski
% R = mychol(A);
% CC = R*R';
% [x3,res3,iter3] = CGprecd(CC,A,b,x0,tol,maxiter);
% err = norm(x1 - xex) ;
% disp(['error Preconditioner CG with Cholewski = ', num2str(err), ', iter = ', num2str(iter3)])

% % Preconditioner CG with Incomplete Cholewski
% SIC = myinccompletechol(A);
% % SIC = ichol(A);
% CIC = SIC*SIC';
% [x4,res4,iter4] = CGprecd(CIC,A,b,x0,tol,maxiter);
% err = norm(x4 - xex) ;
% disp(['error Preconditioner CG with Incomplete Cholewski = ', num2str(err), ', iter = ', num2str(iter4)])

% Preconditioner CG with LU(0) decomposition
W = ilu(A);
M = W*W';
[x51,res51,iter51] = CGprecd(M,A,b,x0,tol,maxiter);
err = norm(x51 - xex) ;
disp(['error Preconditioner CG with L+U-speye(size(A)) = ', num2str(err), ', iter = ', num2str(iter51)])


% Preconditioner CG with LU(0) decomposition
[L,U] = ilu(A);
M = (L*U);
[x5,res5,iter5] = CGprecd(M,A,b,x0,tol,maxiter);
err = norm(x5 - xex) ;
disp(['error Preconditioner CG with LU(0) decomposition = ', num2str(err), ', iter = ', num2str(iter5)])

% Preconditioner CG with LU(1) decomposition
[L,U] = ilu(A*A);
M = (L*U);
[x6,res6,iter6] = CGprecd(M,A,b,x0,tol,maxiter);
err = norm(x6 - xex) ;
disp(['error Preconditioner CG with LU(1) decomposition = ', num2str(err), ', iter = ', num2str(iter6)])




% Visualize convergence versus iterations
s0 = 'Conjugate gradient';
semilogy(res0,'-x', 'DisplayName',s0)
hold on
s1 = 'Preconditioner CG with Gauss Seidel';
semilogy(res1,'-p', 'DisplayName',s1)
hold on
s2 = 'Preconditioner CG with SSOR';
semilogy(res2,'-<', 'DisplayName',s2)
% s3 = 'Preconditioner CG with Cholewski';
% semilogy(res3,'->', 'DisplayName',s3)
% s4 = 'Preconditioner CG with Incomplete Cholewski';
% semilogy(res4,'-*', 'DisplayName',s4)

s51 = 'Preconditioner CG with L+U-speye(size(A))';
semilogy(res51,'-^ ', 'DisplayName',s51)
s5 = 'Preconditioner CG with LU(0) decomposition';
semilogy(res5,'-o ', 'DisplayName',s5)
s6 = 'Preconditioner CG with LU(1) decomposition';
semilogy(res6,'-d ', 'DisplayName',s6)

legend(s0,s1,s2,s51,s5,s6)















