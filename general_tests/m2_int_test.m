clear all
%
%
%
warning('PIM_hankel.m has been changed so that it containsm1 and m2 functions')
%
%
%

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

f = @(x, y) 1i*(besselh(0, k*abs(x - y)) - log(k*abs(x - y)))/4;



for j = 1:length(s)
    
    mat_approx(j) = integral(@(y) f(s(j), y), a, b);

    our_approx(j) = sum(h*m2(k, s(j), t, C1, C2));

%     PIM_fun = PIM_int_hankel_f(k, s(j), h, t, 1, t_grid, C1, C2);

    err1(j) = mat_approx(j) - our_approx(j);

%     err2(j) = mat_approx - PIM_fun;


end

figure()
plot(s, err1, 'DisplayName', 'Modular code', 'LineStyle', '--')
% hold on
% plot(s, err2, 'DisplayName', 'Not modular code', 'LineStyle', '-.')
% legend show
title('Plot of error difference')

figure()
plot(s, mat_approx, 'DisplayName', 'Matlab approximation')
hold on 
plot(s, our_approx, 'DisplayName', 'Our approximation')
hold off
legend show
title('Comparison plot')