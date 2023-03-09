function z = v1_weights_mid(a,b,k,s)
% Computes \int_{a}^{b} ln(k|s - t|) dt.
%
% Problem parameters:
% a is the lower bound of the integral (could be vector)
% b is the upper bound of the integral (must be same size as a)
% k is the wavenumber
% s is a fixed point (scalar)
%
% No need to split - each alternative boils down to the same formula:
z = (s-a).*log(abs(k*(s-a))) + (b-s).*log(abs(k*(b-s))) + a-b;
end