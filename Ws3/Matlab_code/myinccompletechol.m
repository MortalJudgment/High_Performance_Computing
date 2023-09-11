function ChI=myinccompletechol(A)
nn = size(A,1);
ChI = tril(A);
for k=1:nn
    ChI(k,k) = sqrt(ChI(k,k));
    for j=k+1:nn
        if ChI(j,k) ~= 0
            ChI(j,k) = ChI(j,k)/ChI(k,k);
        end
    end
    for j=k+1:nn
        for i=j:nn
            if ChI(i,j) ~= 0
                ChI(i,j) = ChI(i,j)-ChI(i,k)*ChI(j,k);
            end
        end
    end
end