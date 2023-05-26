% script to test that m1 and m2 are working as they should be
%
% We will fix s at a range of values and then comoute a range of values for
% t, between 0 and L and compare to make sure that:
% i*H_0^(1)(k|s - t|)/4 = sigma1(s, t)m1(s, t) + m2(s, t)
clear all

a = 0;
b = 2*pi;
N = 1000;
h = (b - a)/N;

t = [a:h:b];

s = [a - 2, a, a+ 0.01, (a+b)/2, b-0.001, b, b+0.0001, b+20];

err(j) = zeros(size(s));
for j = 1:length(s)
    dist = abs(s(j) - t);

    hankel = 1i*besselh(0, k*dist)/4;
    our_approx = m1(k, s, t, C1, C2).*

    err(j) = sum()

end

