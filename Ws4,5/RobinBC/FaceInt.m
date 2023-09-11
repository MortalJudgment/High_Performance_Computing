function M = FaceInt(M,la,lb,lc,ld,ic1,ic2,nx,ny)

%Reception face
%-------------------------------
lab = (lb-la)/ic1 -1  ;
M.lintr(M.Nintr+1: M.Nintr +lab) = (la+ic1:ic1:lb-ic1) ;

% Outward Normal
M.nxr(M.Nintr+1: M.Nintr +lab) = nx ;
M.nyr(M.Nintr+1: M.Nintr +lab) = ny ;

M.Nintr = M.Nintr + lab;
M.lvoisr(M.Nvois+1) = M.lvoisr(M.Nvois) +lab;


%Send face
%-------------------------------
lcd = (ld-lc)/ic1 -1 ;
M.linte(M.Ninte+1: M.Ninte + lcd) = (lc+ic1:ic1:ld-ic1) ;

% Inward Normal
M.nxe(M.Ninte+1: M.Ninte + lcd) = -nx ;
M.nye(M.Ninte+1: M.Ninte + lcd) = -ny ;

M.Ninte = M.Ninte + lcd;
M.lvoise(M.Nvois+1) = M.lvoise(M.Nvois) + lcd;


