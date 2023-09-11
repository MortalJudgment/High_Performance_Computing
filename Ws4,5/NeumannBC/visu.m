function visu(indv,T,M,P)
if(indv == 1)
    figure(1) ;
    hold on
    hidden off
    for i=1:T.Ndom
        m1 = M(i).Nx ;
        m3 = M(i).Ny ;
        [xx,yy] = meshgrid(M(i).x,M(i).y) ;
        rng = i/T.Ndom ;
        on = ones(m1,m3)*rng;
        mesh(xx',yy',on);  view(2)
    end
    return
end

if(indv == 2)
    figure(2) ;
    hold on
    hidden off
    for i=1:T.Ndom
        m1 = M(i).Nx ;
        m3=M(i).Ny ;
        [xx,yy] = meshgrid(M(i).x,M(i).y) ;
        on = reshape(P(i).uana, m1, m3);        
        surf(xx',yy',on); view(3)
    end
    return
end

if(indv == 3)
    figure(3) ;
    hold on
    hidden off
    for i=1:T.Ndom
        m1=M(i).Nx ;
        m3=M(i).Ny ;
        [xx,yy] = meshgrid(M(i).x,M(i).y) ;
        onn = reshape(P(i).usol, m1, m3);        
        surf(xx',yy',onn); view(3)
    end
    return
end

