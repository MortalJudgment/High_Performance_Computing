function A = createMatrix(N,beta,alpha)  
d=1:1:N;
D = sparse(N,N) ;
for i=1:N
    D(i,i)=alpha*d(i);
end
e = ones(N,1) ;
sb = beta*e ;
S = spdiags( [e, sb], [0  1] , N,N ) ;
 A=S*D*inv(S);