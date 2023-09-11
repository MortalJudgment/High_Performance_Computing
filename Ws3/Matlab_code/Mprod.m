function[y]=Mprod(iter, x,invT,U )  

%M=eye(n,n)+U*(abs(I(m))*inv(T)-eye(rmax,rmax))*U';

if iter == 1
    y = x ;
else
    y = U'*x ;
    y = invT*y -y ;
    y = x+ U*y ;
end
