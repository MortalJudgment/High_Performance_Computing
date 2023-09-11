function Ch = mychol(A)
Ch=tril(A);
nn=size(A,1);
for k=1:nn
    Ch(k,k)=sqrt(Ch(k,k));
    Ch(k+1:nn,k)=Ch(k+1:nn,k)/Ch(k,k);
    for j=k+1:nn
        Ch(j:nn,j)=Ch(j:nn,j)-Ch(j:nn,k)*Ch(j,k);
    end
end