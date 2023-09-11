function [c,s] = rotmat(a,b)
%
% Compute the Givens rotation matrix parameters for a and b.
if ( b == 0.0 ),
      c = 1.0;
      s = 0.0;
else 
      t= sqrt(a^2 + b^2);
      c = a/t;
      s = b/t;
end
