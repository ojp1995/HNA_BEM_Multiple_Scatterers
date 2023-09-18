function z = v1_weights_mid_vec_test(a,b,k,s)
% Computes \int_{a}^{b} ln(k|s - t|) dt.
%
% Problem parameters:
% a is the lower bound of the integral (could be vector)
% b is the upper bound of the integral (must be same size as a)
% k is the wavenumber
% s is a fixed point (scalar)
%
z = zeros(size(a));
select = a~=s & b~=s;  % added to get around the case when a==s and we know it should be zero
% No need to split - each alternative boils down to the same formula:
z(select) = (s-a(select)).*log(abs(k*(s-a(select)))) + (b(select)-s).*log(abs(k*(b(select)-s))) + a(select)-b(select);
end