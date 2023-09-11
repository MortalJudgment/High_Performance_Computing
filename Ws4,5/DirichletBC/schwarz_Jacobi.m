function  [P,res] = schwarz_Jacobi(P,M,T,Mg,epsi,itmax)
vari = zeros(T.Ndom,1) ;

for it=1:itmax
    for i=1:T.Ndom
        P(i).sm = P(i).sm0  ;  % sm0 has the original Dirichlet conditions 
		                       % on the outer bofders
        %--------------------------------------
        %  Updating boundary values
        %--------------------------------------        
        if(M(i).Nintr~= 0)
            for j = 1:M(i).Nvois
                % list of boundary points to be updated
                indSi = (M(i).lvoisr(j)+1 : M(i).lvoisr(j+1));
                indAi = M(i).lintr(indSi) ;
                
 				% List of these points in the neighbour
                idvois = M(i).novois(j) ;
                ij = M(idvois).list(i)  ; % position of these points in the neighbour
                indSj = (M(idvois).lvoise(ij)+1 : M(idvois).lvoise(ij+1)) ;
                indAj = M(idvois).linte(indSj) ;
                P(i).sm(indAi) = P(idvois).usol(indAj);
            end
        end
        %--------------------------------------
        %    End  Updating boundary values
        %--------------------------------------        
        usol0  = P(i).usol ;                %old solution
        P(i).usol = P(i).A\P(i).sm;         %New solution
%         P(i).usol = gmres(P(i).A,P(i).sm);
        vari(i) = max(abs(usol0 - P(i).usol));        
    end
    
    var = max(vari) ;
    res(it) = var;
    sprintf('it= %d  ,var= %e\n', it,var)
%     visu(2,T,M,P) ;
    if( var < epsi ) return;  end
    
end