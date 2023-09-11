function[A] = create_matrix(n,beta,alp)  
d = 1:1:n;
D = sparse(n,n);
for i=1:n
    D(i,i) = alp*d(i);
end

c = ones(n,1) ;
cb = beta*c ;
S = spdiags( [c, cb], [0  1], n, n ) ;

 A = S*D*inv(S);