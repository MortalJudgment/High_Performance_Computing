function [P] = schwarz(P,M,T,Mg,tol,maxit)
it = 1;
error = 1;
err = zeros(1,T.Ndom);
m1 = M.Nx ;
m2 = M.Ny ;
mm = m1*m2;

while it < maxit && error >= tol
    for i = 1:T.Ndom
        P(i).sm = P(i).sm0  ;
        for j = 1:M(i).Ndir         % On Dirichlet boundary condition
            P(i).sm(M(i).ldir(j)) = P(i).udir(M(i).ldir(j));
        end
        
        for j = 1:M(i).Nintr        % Dirichlet interface conditions
            if i < T.Ndom
                P(i).sm(M(i).lintr(j)) = P(i+1).usol(M(i+1).linte(j));
            else
                P(i).sm(M(i).lintr(j)) = P(1).usol(M(i).linte(j));
            end
        end
        oldsol = P(i).usol ;
        P(i).usol = P(i).A\P(i).sm;
        err(i) = norm(oldsol - P(i).usol);
    end
    it = it + 1;
    error = max(err);
end
end


