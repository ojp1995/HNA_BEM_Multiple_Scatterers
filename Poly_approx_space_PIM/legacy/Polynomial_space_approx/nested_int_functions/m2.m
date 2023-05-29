function z = m2(k, s, C1, C2)
% This function is the smooth part of the function from splitting k(s) =
% sigmna1(s)m1(s) + sigma2(s)m2(s), with the smoothing function applied.
% Here m2(s) = k(s) - sigma1(s)m1_tilde(s)
%
% Problem parameters:
% k is the wavenumber
% s is the input of the function
% C1 is the first constant for the smoothing funciton, when k*s<C1,
% smoothing function = 1.
% C2 is the second constant for the smoothing function, when k*s>= C2,
% smoothing function = 0.
% 
% Discretisation parameters:

% sigma1(s) = ln(k|s|);

z = 1i*besselh(0, k*s)/4 - m1_tilde(k, s, C1, C2).*log(k*s);