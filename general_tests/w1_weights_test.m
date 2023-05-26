% in this function we will be testing that the integrator w1_weights.m
% works, by computing the integral \int_a^b ln(k|s - t|) dt, for a range of
% values of s. To compare this we will be using the matlab inbuilt
% integrator.

clear all

addpath('/Users/ojp18/Dropbox/Mac/Documents/GitHub/HNA_BEM_Multiple_Scatterers/General_functions')

k = 0.85;
lambda = 2*pi/k;
a = 1.35;
b = 1 + lambda/10;

s = linspace(0, 10, 1000);

N = 2^20;
h = (b - a)/N;
C1 = 1;
C2 = 2*pi;

t = [a+h/2:h:b- h/2];
t_grid = [a:h:b];


err = zeros(size(s));

f = @(x, y) log(k*abs(x - y));

for j = 1:length(s)

    matlab_val(j) = integral(@(y) f(s(j), y), a, b);

    our_approx(j) = sum(w1_weights(k, s(j), t_grid(1:end - 1), t_grid(2:end)));

    err(j) = sum(matlab_val(j) - our_approx(j));

end

% err
% 
% [matlab_val.' our_approx.']

figure()
plot(s, err)
title('Plot of the difference between Matlabs integrator and our approximation')

