function z = m1_tilde( k, s, C1, C2)
% This function computes m1(s) = \chi(ks)*(-1/2*pi)*J_{0}(ks), where J_{0}
% is the bessel function of the first kind.
%
% Problem parameters:
% k is the wavenmuber
% s is the input of the function
% C1 is the first constant for the smoothing funciton, when k*s<C1,
% smoothing function = 1.
% C2 is the second constant for the smoothing function, when k*s>= C2,
% smoothing function = 0.

 z = -smoothing_fun(k*s, C1, C2).*besselj(0, k*s)/(2*pi);
 
end
