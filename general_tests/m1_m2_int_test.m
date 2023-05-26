% m1_m2_int test
% 

clear all

addpath('/Users/ojp18/Dropbox/Mac/Documents/GitHub/HNA_BEM_Multiple_Scatterers/General_functions')
k = 0.85;
lambda = 2*pi/k;
a = 1.35;
b = 1 + lambda/10;

s = linspace(0, 10, 1000);

N = 2^10;
h = (b - a)/N;
C1 = 1;
C2 = 2*pi;

t = [a+h/2: h: b - h/2];
t_grid = [a:h:b];

err1 = zeros(size(s));
err2 = zeros(size(s));

f = @(x, y) 1i*besselh(0, k*abs(x - y))/4;

for j = 1:length(s)
    
    mat_approx = integral(@(y) f(s(j), y), a, b);

    our_approx = sum(m1(k, s(j), t, C1, C2).*w1_weights(k, s(j), ...
        t_grid(1:end - 1), t_grid(2:end))...
        + h*m2(k, s(j), t, C1, C2));

    PIM_fun = PIM_int_hankel_f(k, s(j), h, t, 1, t_grid, C1, C2);

    err1(j) = mat_approx - our_approx;

    err2(j) = mat_approx - PIM_fun;


end

figure()
plot(s, err1, 'DisplayName', 'Modular code', 'LineStyle', '--')
hold on
plot(s, err2, 'DisplayName', 'Not modular code', 'LineStyle', '-.')
legend show
title('Plot of error difference')